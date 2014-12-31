module GSLIntegrationPort where

import Data.List as L

minPositiveFloat :: RealFloat a => a -> a
minPositiveFloat a = encodeFloat 1 $ fst (floatRange a) - floatDigits a

machEsp :: RealFrac a => a
machEsp = until p (/2) 1
  where p eps = eps/2 + 1 == 1

dbl_epsilon = machEsp :: Double
dbl_min     = minPositiveFloat 0.0 :: Double

abl :: [a] -> [a]
abl = reverse . tail . reverse

dotprod :: [Double] -> [Double] -> Double
dotprod x y  = sum $ zipWith (*) x y

sjb_qng :: (Double -> Double) -> Double -> Double -> (Double, Double)
sjb_qng f a b = (res, err)
  where
    half_length     =  0.5 * (b - a)
    abs_half_length = abs half_length
    center          = 0.5 * (b + a)
    f_center        =  f center
    f_center_abs    = abs f_center
    
    abscissa_x1 = map  (* half_length) x1
    abscissa_x2 = map  (* half_length) x2
    abscissa_x3 = map  (* half_length) x3
    abscissa_x4 = map  (* half_length) x4


    fv1       = map (\a -> f (center + a) ) abscissa_x1
    fv2       = map (\a -> f (center - a) ) abscissa_x1    
    fval_x1   = zipWith (+) fv1 fv2 :: [Double]

    fv3       = map (\a -> f (center + a) ) abscissa_x2
    fv4       = map (\a -> f (center - a) ) abscissa_x2
    fval_x2   = zipWith (+) fv3 fv4 :: [Double]

    fv5       = map (\a -> f (center + a) ) abscissa_x3
    fv6       = map (\a -> f (center - a) ) abscissa_x3
    fval_x3   = zipWith (+) fv5 fv6 :: [Double]

    fv7       = map (\a -> f (center + a) ) abscissa_x4
    fv8       = map (\a -> f (center - a) ) abscissa_x4
    fval_x4   = zipWith (+) fv7 fv8 :: [Double]

    
    -- res10
    res10    = dotprod fval_x1 w10

    -- res21
    res21_fc = (last w21b) * f_center
    res21_x1 = dotprod fval_x1 w21a
    res21_x2 = dotprod fval_x2 w21b
    res21    = res21_fc + res21_x1 + res21_x2

    -- res43
    res43_fc = (last w43b) * f_center
    res43_a  = dotprod (fval_x1 ++ fval_x2) w43a
    res43_b  = dotprod fval_x3 w43b
    res43    = res43_fc + res43_a + res43_b

    -- res87
    res87_fc = (last w87b) * f_center
    res87_a  = dotprod (fval_x1 ++ fval_x2 ++ fval_x3) w87a
    res87_b  = dotprod fval_x4 w87b
    res87    = res87_fc + res87_a + res87_b

    res = res87 * half_length

    
    -- resabs
    resabs_fc = (last w21b) * f_center
    resabs_x1 = wtabs fv1 fv2 w21a
    resabs_x2 = wtabs fv3 fv4 w21b
    resabs    = abs_half_length * (resabs_fc + resabs_x1 + resabs_x2)

    -- resasc
    mean        = 0.5 * res21
    resasc_fc   = (last w21b) * (abs (f_center - mean))
    resasc_main = sum $ zipWith6 (\a b c d w1 w2 -> (w1 * ((abs (a - mean)) + (abs (b - mean)))) + (w2 * ((abs (c - mean)) + (abs (d - mean))))) fv1 fv2 fv3 fv4 w21a w21b
    resasc      = (resasc_fc + resasc_main) * abs_half_length

    -- result_konrod
    result_konrod = res87 * half_length

    err = rescale_error ((res87 - res43) * half_length) resabs resasc


wtabs al bl wtl = sum $ zipWith3 (\a b wt -> ((abs a) + (abs b)) * wt ) al bl wtl

rescale_error :: Double -> Double -> Double -> Double
rescale_error err_in result_abs result_asc  = err_out
  where
    err   = abs err_in
    base  = (200.0 * err) / result_asc
    scale = base ** 1.5
    min_err = 50.0 * dbl_epsilon * result_abs

    err2  = if ( (result_asc /= 0)  && (err /= 0)) then
              if (scale < 1) then (result_asc * scale) else result_asc
            else err

    err_out = if (result_abs > dbl_min / (50 * dbl_epsilon)) then
                if (min_err > err) then min_err else err2
              else err2

-- abcissa weights

x1 = [  0.973906528517171720077964012084452,
        0.865063366688984510732096688423493,
        0.679409568299024406234327365114874,
        0.433395394129247190799265943165784,
        0.148874338981631210884826001129720
     ]


-- /* w10, weights of the 10-point formula */
w10 = [  0.066671344308688137593568809893332,
         0.149451349150580593145776339657697,
         0.219086362515982043995534934228163,
         0.269266719309996355091226921569469,
         0.295524224714752870173892994651338
      ]

-- /* x2, abscissae common to the 21-, 43- and 87-point rule */
x2 = [  0.995657163025808080735527280689003,
        0.930157491355708226001207180059508,
        0.780817726586416897063717578345042,
        0.562757134668604683339000099272694,
        0.294392862701460198131126603103866
     ]

-- /* w21a, weights of the 21-point formula for abscissae x1 */
w21a = [  0.032558162307964727478818972459390,
          0.075039674810919952767043140916190,
          0.109387158802297641899210590325805,
          0.134709217311473325928054001771707,
          0.147739104901338491374841515972068
       ]

-- /* w21b, weights of the 21-point formula for abscissae x2 */
w21b = [  0.011694638867371874278064396062192,
          0.054755896574351996031381300244580,
          0.093125454583697605535065465083366,
          0.123491976262065851077958109831074,
          0.142775938577060080797094273138717,
          0.149445554002916905664936468389821
       ]

-- /* x3, abscissae common to the 43- and 87-point rule */
x3 = [  0.999333360901932081394099323919911,
        0.987433402908088869795961478381209,
        0.954807934814266299257919200290473,
        0.900148695748328293625099494069092,
        0.825198314983114150847066732588520,
        0.732148388989304982612354848755461,
        0.622847970537725238641159120344323,
        0.499479574071056499952214885499755,
        0.364901661346580768043989548502644,
        0.222254919776601296498260928066212,
        0.074650617461383322043914435796506
     ]

-- /* w43a, weights of the 43-point formula for abscissae x1, x3 */
w43a = [  0.016296734289666564924281974617663,
          0.037522876120869501461613795898115,
          0.054694902058255442147212685465005,
          0.067355414609478086075553166302174,
          0.073870199632393953432140695251367,
          0.005768556059769796184184327908655,
          0.027371890593248842081276069289151,
          0.046560826910428830743339154433824,
          0.061744995201442564496240336030883,
          0.071387267268693397768559114425516
       ]

-- /* w43b, weights of the 43-point formula for abscissae x3 */
w43b = [  0.001844477640212414100389106552965,
          0.010798689585891651740465406741293,
          0.021895363867795428102523123075149,
          0.032597463975345689443882222526137,
          0.042163137935191811847627924327955,
          0.050741939600184577780189020092084,
          0.058379395542619248375475369330206,
          0.064746404951445885544689259517511,
          0.069566197912356484528633315038405,
          0.072824441471833208150939535192842,
          0.074507751014175118273571813842889,
          0.074722147517403005594425168280423
       ]

-- /* x4, abscissae of the 87-point rule */
x4 = [  0.999902977262729234490529830591582,
        0.997989895986678745427496322365960,
        0.992175497860687222808523352251425,
        0.981358163572712773571916941623894,
        0.965057623858384619128284110607926,
        0.943167613133670596816416634507426,
        0.915806414685507209591826430720050,
        0.883221657771316501372117548744163,
        0.845710748462415666605902011504855,
        0.803557658035230982788739474980964,
        0.757005730685495558328942793432020,
        0.706273209787321819824094274740840,
        0.651589466501177922534422205016736,
        0.593223374057961088875273770349144,
        0.531493605970831932285268948562671,
        0.466763623042022844871966781659270,
        0.399424847859218804732101665817923,
        0.329874877106188288265053371824597,
        0.258503559202161551802280975429025,
        0.185695396568346652015917141167606,
        0.111842213179907468172398359241362,
        0.037352123394619870814998165437704
     ]

-- /* w87a, weights of the 87-point formula for abscissae x1, x2, x3 */
w87a = [  0.008148377384149172900002878448190,
          0.018761438201562822243935059003794,
          0.027347451050052286161582829741283,
          0.033677707311637930046581056957588,
          0.036935099820427907614589586742499,
          0.002884872430211530501334156248695,
          0.013685946022712701888950035273128,
          0.023280413502888311123409291030404,
          0.030872497611713358675466394126442,
          0.035693633639418770719351355457044,
          0.000915283345202241360843392549948,
          0.005399280219300471367738743391053,
          0.010947679601118931134327826856808,
          0.016298731696787335262665703223280,
          0.021081568889203835112433060188190,
          0.025370969769253827243467999831710,
          0.029189697756475752501446154084920,
          0.032373202467202789685788194889595,
          0.034783098950365142750781997949596,
          0.036412220731351787562801163687577,
          0.037253875503047708539592001191226
       ]

-- /* w87b, weights of the 87-point formula for abscissae x4    */
w87b = [  0.000274145563762072350016527092881,
          0.001807124155057942948341311753254,
          0.004096869282759164864458070683480,
          0.006758290051847378699816577897424,
          0.009549957672201646536053581325377,
          0.012329447652244853694626639963780,
          0.015010447346388952376697286041943,
          0.017548967986243191099665352925900,
          0.019938037786440888202278192730714,
          0.022194935961012286796332102959499,
          0.024339147126000805470360647041454,
          0.026374505414839207241503786552615,
          0.028286910788771200659968002987960,
          0.030052581128092695322521110347341,
          0.031646751371439929404586051078883,
          0.033050413419978503290785944862689,
          0.034255099704226061787082821046821,
          0.035262412660156681033782717998428,
          0.036076989622888701185500318003895,
          0.036698604498456094498018047441094,
          0.037120549269832576114119958413599,
          0.037334228751935040321235449094698,
          0.037361073762679023410321241766599
       ]

--------------------------------------------------------------------------------



data GKStruct = GKStruct {
   result :: Double
 , abserr :: Double
 , resabs :: Double
 , resasc :: Double  
}

integration_qk :: 
     Int        -- ^ n 
 ->  [Double]   -- ^ xgk 
 ->  [Double]   -- ^ wg 
 ->  [Double]   -- ^ wgk 
 ->  [Double]   -- ^ fv1 
 ->  [Double]   -- ^ fv2
 ->  (Double -> Double) -- ^ f
 ->  Double     -- ^ a
 ->  Double     -- ^ b
 ->  GKStruct
integration_qk n xgk wg wgk fv1 fv2 f a b = gsk
  where
    center          = 0.5 * (a + b)
    half_length     = 0.5 * (b - a)
    abs_half_length = abs (half_length)
    f_center        = f center

    -- NOTE: when idx is used, it is only when 
    nd              = fromIntegral n :: Double
    idx             = round $ (nd / 2.0) - 1.0  
    jmax_1          = (nd - 1.0) / 2.0
    jmax_2          = nd / 2.0

    -- jvec_1 = takeWhile ( < jmax_1) [0..] :: [Int]
    -- jvec_2 = takeWhile ( < jmax_2) [0..] :: [Int]

    jvec_1     = takeWhile ( <= (floor jmax_1)) [0..] :: [Int]
    jvec_2     = takeWhile ( <= (floor jmax_2)) [0..] :: [Int]
    jtw_vec    = map (\j -> j * 2 + 1) jvec_1
    jtwm1_vec  = map (\j -> j * 2)     jvec_2

    abcissa_j1 = map (\i -> half_length * (xgk !! i )) jtw_vec
    abcissa_j2 = map (\i -> half_length * (xgk !! i )) jtwm1_vec
    fv1_j1     = map (\x -> f (center - x)) abcissa_j1
    fv2_j1     = map (\x -> f (center + x)) abcissa_j1
    fv1_j2     = map (\x -> f (center - x)) abcissa_j2
    fv2_j2     = map (\x -> f (center + x)) abcissa_j2    
    fsum_j1    = zipWith (+) fv1_j1 fv2_j1
    fsum_j2    = zipWith (+) fv1_j2 fv2_j2

    -- Gauss result
    result_gauss_n  = if (mod n 2) == 0 then (f_center * (wg !! idx)) else 0.0
    result_gauss_j1 = dotprod (map (wg !! ) jvec_1) fsum_j1
    result_gauss    = result_gauss_n + result_gauss_j1

    -- Konrod result
    konrod_fc = f_center * (wgk !! (n-1))
    konrod_j1 = dotprod (map (wgk !!) jtw_vec)   fsum_j1
    konrod_j2 = dotprod (map (wgk !!) jtwm1_vec) fsum_j2
    konrod    = konrod_fc + konrod_j1 + konrod_j2
    result_konrod = konrod * half_length

    -- Abs Result
    abs_fc = abs konrod_fc
    abs_j1 = dotprod (map (wgk !!) jtw_vec)   (abssum fv1_j1 fv2_j1)
    abs_j2 = dotprod (map (wgk !!) jtwm1_vec) (abssum fv1_j2 fv2_j2)
    result_abs    = abs_half_length * (abs_fc + abs_j1 + abs_j2)
    
    -- Asc result
    fv1_abs_demean = map abs $ map (flip (-) center) $ weave fv1_j2 fv1_j1
    fv2_abs_demean = map abs $ map (flip (-) center) $ weave fv2_j2 fv2_j1
    mean       = 0.5 * result_konrod
    asc_fc     = (wgk !! (n-1)) * (abs (f_center - mean))
    asc_1      = dotprod (take (n - 2) wgk) $ zipWith (+) fv1_abs_demean fv2_abs_demean
    result_asc = abs_half_length * (asc_fc + asc_1)

    -- Abs err
    err    = (result_konrod - result_gauss) * half_length
    abserr = rescale_error err result_abs result_asc

    -- Construct error
    gsk    = GKStruct result_konrod abserr result_abs result_asc
    

weave :: [a] -> [a] -> [a]
weave x []          = x
weave [] y          = y
weave (x:xs) (y:ys) = x:y:(weave xs ys)

abssum :: [Double] -> [Double] -> [Double]
abssum xlist ylist = zipWith (\x y -> (abs x) + (abs y)) xlist ylist


-- qags integration
-- port of GSL's facillities

-- qags ::    (Double -> Double)  -- ^ integrand
--         -> Double              -- ^ lower bound
--         -> Double              -- ^ upper bound
--         -> Double              -- ^ epsabs
--         -> Double              -- ^ epsrel
--         -> Double              -- ^ limit
--         -> (Double -> Double)  -- ^ gsl_integration_rule
--         -> (Double, Double)    -- ^ (result, error)
-- qags f a b epsabs epsrel limit gsl_integration_rule = (result, error)
--   where
    
