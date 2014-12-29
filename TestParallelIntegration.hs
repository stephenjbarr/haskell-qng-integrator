import Control.Seq as Seq
import Control.Monad.Par
import Control.DeepSeq
import GSLIntegrationPort
import Data.List
import qualified Numeric.GSL.Integration as GSL
import Criterion.Main
  
parMapChunk :: (Eq b, NFData b) => Int -> (a -> b) -> [a] -> Par [b]
parMapChunk n f xs = fmap concat $ parMap (map f) (chunk n xs)

chunk :: Int -> [a] -> [[a]]
chunk _ [] = []
chunk n xs = as : chunk n bs where (as,bs) = splitAt n xs

ncores = 12

quad_finite = GSL.integrateQNG   1E-6

cdiv q d = ceiling $ (fromIntegral q) / (fromIntegral d)

avg l = let (t,n) = foldl' (\(b,c) a -> (a+b,c+1)) (0,0) l 
        in realToFrac(t)/realToFrac(n)

normpdf :: Double -> Double -> Double -> Double
normpdf x mu sigma = y  * myexp where
  y = (1.0/(sigma*(2.0*pi)**(0.5)))
  myexp = exp( -1 * (x - mu)**(2.0) / (2.0 * sigma**2.0))

stdnorm :: Double -> Double
stdnorm x = normpdf x 0.0 1.0



-- linspace :: (Num a, Integral b)  => a  -> a -> b -> [a]
linspace :: Double -> Double -> Int -> [Double]
linspace lb ub npts = pts
  where
    nxi   = npts - 1
    nx    = fromIntegral nxi
    width = (ub - lb) / nx 
    pts   = map (\i -> lb + (fromIntegral i) * width) [0..nxi]
    

-- moving region of width 1 from (-3,-2) to (2,3)
lb_ub_list :: Int -> [(Double, Double)]
lb_ub_list ni = lbub
  where
    lbs  = linspace (-3.0) 2.0 ni
    ubs  = linspace (-2.0) 3.0 ni
    lbub = zip lbs ubs


-- run test
sequential_test_sjb :: Int -> Double    
sequential_test_sjb n = avg intlist
  where
    intf (a,b)  = fst $ sjb_qng stdnorm a b
    intlist     = map intf (lb_ub_list n)


sequential_test_gsl :: Int -> Double
sequential_test_gsl n = avg intlist
  where
    intf (a,b)  = fst $ quad_finite stdnorm a b
    intlist     = map intf (lb_ub_list n)


parallel_test_gsl :: Int -> Double
parallel_test_gsl n = avg intlist
  where
    intf (a,b)  = fst $ quad_finite stdnorm a b
    intlist     = runPar $  parMap intf (lb_ub_list n)

parallel_test_gsl_chunk :: Int -> Double
parallel_test_gsl_chunk n = avg intlist
  where
    chunk_size = cdiv n ncores
    intf (a,b)  = fst $ quad_finite stdnorm a b
    intlist     = runPar $  parMapChunk chunk_size intf (lb_ub_list n)
    


parallel_test_sjb :: Int -> Double
parallel_test_sjb n = avg intlist
  where
    intf (a,b)  = fst $ sjb_qng stdnorm a b
    intlist     = runPar $  parMap intf (lb_ub_list n)

parallel_test_sjb_chunk :: Int -> Double
parallel_test_sjb_chunk n = avg intlist
  where
    chunk_size = cdiv n ncores
    intf (a,b)  = fst $ sjb_qng stdnorm a b
    intlist     = runPar $  parMapChunk chunk_size intf (lb_ub_list n)



main = defaultMain [

      bgroup "sequential_sjb" [ bench "10"    $ whnf sequential_test_sjb 10
                              , bench "100"   $ whnf sequential_test_sjb 100
                              , bench "1000"  $ whnf sequential_test_sjb 1000
                              , bench "10000" $ whnf sequential_test_sjb 10000
                              ]
      ,

      bgroup "sequential_gsl" [ bench "10"    $ whnf sequential_test_gsl 10     
                              , bench "100"   $ whnf sequential_test_gsl 100   
                              , bench "1000"  $ whnf sequential_test_gsl 1000 
                              , bench "10000" $ whnf sequential_test_gsl 10000
                              ]
      ,

      bgroup "parallel_sjb" [ bench "10"    $ whnf parallel_test_sjb (ncores * 10    )
                            , bench "100"   $ whnf parallel_test_sjb (ncores * 100   )
                            , bench "1000"  $ whnf parallel_test_sjb (ncores * 1000  )
                            , bench "10000" $ whnf parallel_test_sjb (ncores * 10000 )
                            ]
      ,

      bgroup "parallel_gsl" [ bench "10"    $ whnf parallel_test_gsl   (ncores * 10    ) 
                            , bench "100"   $ whnf parallel_test_gsl (ncores * 100   )
                            , bench "1000"  $ whnf parallel_test_gsl (ncores * 1000  )
                            , bench "10000" $ whnf parallel_test_gsl (ncores * 10000 )
                            ]

      ,

      bgroup "parallel_sjb_chunk" [ bench "10"    $ whnf parallel_test_sjb_chunk (ncores * 10    )
                                  , bench "100"   $ whnf parallel_test_sjb_chunk (ncores * 100   )
                                  , bench "1000"  $ whnf parallel_test_sjb_chunk (ncores * 1000  )
                                  , bench "10000" $ whnf parallel_test_sjb_chunk (ncores * 10000 )
                                  ]
      ,

      bgroup "parallel_gsl_chunk" [ bench "10"    $ whnf parallel_test_gsl_chunk   (ncores * 10    ) 
                                  , bench "100"   $ whnf parallel_test_gsl_chunk (ncores * 100   )
                                  , bench "1000"  $ whnf parallel_test_gsl_chunk (ncores * 1000  )
                                  , bench "10000" $ whnf parallel_test_gsl_chunk (ncores * 10000 )
                                  ]


      ]
