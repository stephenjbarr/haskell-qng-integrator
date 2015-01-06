import Control.Seq as Seq
import Control.Monad.Par
import Control.DeepSeq
import GSLIntegrationPort
import Data.List
import qualified Numeric.GSL.Integration as GSL
import Criterion.Main
import qualified Numeric.Integration.TanhSinh as TS  

parMapChunk :: (Eq b, NFData b) => Int -> (a -> b) -> [a] -> Par [b]
parMapChunk n f xs = fmap concat $ parMap (map f) (chunk n xs)

chunk :: Int -> [a] -> [[a]]
chunk _ [] = []
chunk n xs = as : chunk n bs where (as,bs) = splitAt n xs

ncores = 12



cdiv :: Int -> Int -> Int
cdiv q d = ceiling $ (fromIntegral q) / (fromIntegral d)

avg l = let (t,n) = foldl' (\(b,c) a -> (a+b,c+1)) (0,0) l 
        in realToFrac(t)/realToFrac(n)

normpdf :: Double -> Double -> Double -> Double
normpdf x mu sigma = y  * myexp where
  y = (1.0/(sigma*(2.0*pi)**(0.5)))
  myexp = exp( -1 * (x - mu)**(2.0) / (2.0 * sigma**2.0))

stdnorm :: Double -> Double
stdnorm x = normpdf x 0.0 1.0


test_tol = 1E-6


-- linspace :: (Num a, Integral b)  => a  -> a -> b -> [a]
linspace :: Double -> Double -> Int -> [Double]
linspace lb ub npts = pts
  where
    nxi   = npts - 1
    nx    = fromIntegral nxi
    width = (ub - lb) / nx 
    pts   = map (\i -> lb + (fromIntegral i) * width) [0..nxi]
    

-- moving region of width 1 from (-3,-2) to (2,3)
-- lb_ub_list :: Int -> [(Double, Double)]
-- lb_ub_list ni = lbub
--   where
--     lbs  = linspace (-3.0) 2.0 ni
--     ubs  = linspace (-2.0) 3.0 ni
--     lbub = zip lbs ubs


lb_list :: Int -> [Double]
lb_list n = linspace (-3.0) 2.0 n

rtI :: (Double -> Double) -> Double -> Double
rtI f lb = TS.result $  TS.absolute test_tol $  rtTail TS.trap stdnorm lb

quad_rt     = GSL.integrateQAGIU test_tol 1000
quad_finite = GSL.integrateQNG   test_tol


-- -- run test
sequential_test_trap :: Int -> Double    
sequential_test_trap n = avg intlist
  where
    intlist     = map (rtI stdnorm) (lb_list n)


-- -- run test
sequential_test_qagiu :: Int -> Double    
sequential_test_qagiu n = avg intlist
  where
    intlist     = map (\x -> (fst (quad_rt stdnorm x))) (lb_list n)



parallel_test_trap :: Int -> Double
parallel_test_trap n = avg intlist
  where
    intlist     = runPar $ parMap (rtI stdnorm) (lb_list n)


parallel_test_qagiu :: Int -> Double    
parallel_test_qagiu n = avg intlist
  where
    intlist     = runPar $ parMap (\x -> (fst (quad_rt stdnorm x))) (lb_list n)

parallel_test_trap_chunk :: Int -> Double
parallel_test_trap_chunk n = avg intlist
  where
    chunk_size  = cdiv n ncores
    intlist     = runPar $  parMapChunk chunk_size (rtI stdnorm) (lb_list n)


parallel_test_qagiu_chunk :: Int -> Double
parallel_test_qagiu_chunk n = avg intlist
  where
    chunk_size  = cdiv n ncores
    intlist     = runPar $  parMapChunk chunk_size (\x -> (fst (quad_rt stdnorm x))) (lb_list n)



-- | Integrate from lb to infinity using the x = t/(1-t) substitution
rtTail :: ((Double -> Double) -> Double -> Double -> r) -> (Double -> Double) -> Double -> r
rtTail method f lb  = method (\t -> f(t/(1-t))/square(1-t)) (vsub_inv lb) 1
  where
    square x = x * x

-- | Integrate from (- infinity) to ub using the x = t/(1-t) substitution
-- lftTail :: ((Double -> Double) -> Double -> Double -> r) -> (Double -> Double) -> Double -> r
-- lftTail method f ub  = method (\t -> f(t/(1-t))/square(1-t)) (vsub_inv ub) 1.0 
--   where
--     square x = x * x



poscomp   x = max x 0.0
vsub t      = t/(1-t)
vsub_inv  x = x / (x+1)
square    x = x * x


main = defaultMain [

      bgroup "sequential_trap" [ bench "10"    $ whnf sequential_test_trap 10
                              , bench "100"   $ whnf sequential_test_trap 100
                              , bench "1000"  $ whnf sequential_test_trap 1000
                              , bench "10000" $ whnf sequential_test_trap 10000
                              ]
      ,

      bgroup "sequential_qagiu" [ bench "10"    $ whnf sequential_test_qagiu 10     
                              , bench "100"   $ whnf sequential_test_qagiu 100   
                              , bench "1000"  $ whnf sequential_test_qagiu 1000 
                              , bench "10000" $ whnf sequential_test_qagiu 10000
                              ]
      ,

      bgroup "parallel_trap" [ bench "10"    $ whnf parallel_test_trap (ncores * 10    )
                            , bench "100"   $ whnf parallel_test_trap (ncores * 100   )
                            , bench "1000"  $ whnf parallel_test_trap (ncores * 1000  )
                            , bench "10000" $ whnf parallel_test_trap (ncores * 10000 )
                            ]
      ,

      bgroup "parallel_qagiu" [ bench "10"    $ whnf parallel_test_qagiu   (ncores * 10    ) 
                            , bench "100"   $ whnf parallel_test_qagiu (ncores * 100   )
                            , bench "1000"  $ whnf parallel_test_qagiu (ncores * 1000  )
                            , bench "10000" $ whnf parallel_test_qagiu (ncores * 10000 )
                            ]

      ,

      bgroup "parallel_trap_chunk" [ bench "10"    $ whnf parallel_test_trap_chunk (ncores * 10    )
                                  , bench "100"   $ whnf parallel_test_trap_chunk (ncores * 100   )
                                  , bench "1000"  $ whnf parallel_test_trap_chunk (ncores * 1000  )
                                  , bench "10000" $ whnf parallel_test_trap_chunk (ncores * 10000 )
                                  ]
      ,

      bgroup "parallel_qagiu_chunk" [ bench "10"    $ whnf parallel_test_qagiu_chunk   (ncores * 10    ) 
                                  , bench "100"   $ whnf parallel_test_qagiu_chunk (ncores * 100   )
                                  , bench "1000"  $ whnf parallel_test_qagiu_chunk (ncores * 1000  )
                                  , bench "10000" $ whnf parallel_test_qagiu_chunk (ncores * 10000 )
                                  ]


      ]
