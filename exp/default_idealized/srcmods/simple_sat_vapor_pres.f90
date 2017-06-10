
module simple_sat_vapor_pres_mod

!-----------------------------------------------------------------------
!
!                 simple saturation vapor pressure
!
!      routines for computing the saturation vapor pressure (es) and
!      the derivative of es with respect to temperature
!
!      uses a simple formula based on constant latent heat
!
!-----------------------------------------------------------------------
!
!                               usage
!                               -----
!
!                       call lookup_es  (temp,es)
!
!                       call lookup_des (temp,des)
!
!    arguments
!    ---------
!      temp    intent in       temperature in degrees kelvin
!      es      intent out      saturation vapor pressure in Pascals
!      des     intent out      derivative of saturation vapor pressure
!                              with respect to temperature (Pascals/degree)
!
!-----------------------------------------------------------------------

 use        fms_mod, only:  write_version_number,   &
                            error_mesg, FATAL
 !
 ! use      constants_mod, only:  hlv,rvgas


implicit none
private

 public :: lookup_es, lookup_des
 public :: escomp, descomp ! for backward compatibility
                           ! use lookup_es, lookup_des instead



! ! !-----------------------------------------------------------------------
 interface lookup_es
   module procedure lookup_es_0d, lookup_es_1d, lookup_es_2d, lookup_es_3d
 end interface
 ! for backward compatibility (to be removed soon)
 interface escomp
  module procedure lookup_es_0d, lookup_es_1d, lookup_es_2d, lookup_es_3d
end interface
!-----------------------------------------------------------------------
interface lookup_des
  module procedure lookup_des_0d, lookup_des_1d, lookup_des_2d, lookup_des_3d
end interface
! for backward compatibility (to be removed soon)
interface descomp
  module procedure lookup_des_0d, lookup_des_1d, lookup_des_2d, lookup_des_3d
end interface
! !-----------------------------------------------------------------------
! !  cvs version and tag name
!
! character(len=128) :: version = '$Id: simple_sat_vapor_pres.f90,v 1.6 2002/07/16 22:56:58 fms Exp $'
! character(len=128) :: tagname = '$Name: havana $'
!
! !-----------------------------------------------------------------------
! module variables
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Pre-computed saturation vapor pressure from PyCLES lookup table

real, parameter :: t_table(512) = (/100.15, 100.69765166, 101.24530333, 101.79295499, 102.34060665,&
102.88825832, 103.43590998, 103.98356164, 104.53121331, 105.07886497,&
105.62651663, 106.1741683,  106.72181996, 107.26947162, 107.81712329,&
108.36477495, 108.91242661, 109.46007828, 110.00772994, 110.5553816,&
111.10303327, 111.65068493, 112.19833659, 112.74598826, 113.29363992,&
113.84129159, 114.38894325, 114.93659491, 115.48424658, 116.03189824,&
116.5795499,  117.12720157, 117.67485323, 118.22250489, 118.77015656,&
119.31780822, 119.86545988, 120.41311155, 120.96076321, 121.50841487,&
122.05606654, 122.6037182,  123.15136986, 123.69902153, 124.24667319,&
124.79432485, 125.34197652, 125.88962818, 126.43727984, 126.98493151,&
127.53258317, 128.08023483, 128.6278865,  129.17553816, 129.72318982,&
130.27084149, 130.81849315, 131.36614481, 131.91379648, 132.46144814,&
133.0090998,  133.55675147, 134.10440313, 134.65205479, 135.19970646,&
135.74735812, 136.29500978, 136.84266145, 137.39031311, 137.93796477,&
138.48561644, 139.0332681,  139.58091977, 140.12857143, 140.67622309,&
141.22387476, 141.77152642, 142.31917808, 142.86682975, 143.41448141,&
143.96213307, 144.50978474, 145.0574364,  145.60508806, 146.15273973,&
146.70039139, 147.24804305, 147.79569472, 148.34334638, 148.89099804,&
149.43864971, 149.98630137, 150.53395303, 151.0816047,  151.62925636,&
152.17690802, 152.72455969, 153.27221135, 153.81986301, 154.36751468,&
154.91516634, 155.462818,   156.01046967, 156.55812133, 157.10577299,&
157.65342466, 158.20107632, 158.74872798, 159.29637965, 159.84403131,&
160.39168297, 160.93933464, 161.4869863,  162.03463796, 162.58228963,&
163.12994129, 163.67759295, 164.22524462, 164.77289628, 165.32054795,&
165.86819961, 166.41585127, 166.96350294, 167.5111546,  168.05880626,&
168.60645793, 169.15410959, 169.70176125, 170.24941292, 170.79706458,&
171.34471624, 171.89236791, 172.44001957, 172.98767123, 173.5353229,&
174.08297456, 174.63062622, 175.17827789, 175.72592955, 176.27358121,&
176.82123288, 177.36888454, 177.9165362,  178.46418787, 179.01183953,&
179.55949119, 180.10714286, 180.65479452, 181.20244618, 181.75009785,&
182.29774951, 182.84540117, 183.39305284, 183.9407045,  184.48835616,&
185.03600783, 185.58365949, 186.13131115, 186.67896282, 187.22661448,&
187.77426614, 188.32191781, 188.86956947, 189.41722114, 189.9648728,&
190.51252446, 191.06017613, 191.60782779, 192.15547945, 192.70313112,&
193.25078278, 193.79843444, 194.34608611, 194.89373777, 195.44138943,&
195.9890411,  196.53669276, 197.08434442, 197.63199609, 198.17964775,&
198.72729941, 199.27495108, 199.82260274, 200.3702544,  200.91790607,&
201.46555773, 202.01320939, 202.56086106, 203.10851272, 203.65616438,&
204.20381605, 204.75146771, 205.29911937, 205.84677104, 206.3944227,&
206.94207436, 207.48972603, 208.03737769, 208.58502935, 209.13268102,&
209.68033268, 210.22798434, 210.77563601, 211.32328767, 211.87093933,&
212.418591,   212.96624266, 213.51389432, 214.06154599, 214.60919765,&
215.15684932, 215.70450098, 216.25215264, 216.79980431, 217.34745597,&
217.89510763, 218.4427593,  218.99041096, 219.53806262, 220.08571429,&
220.63336595, 221.18101761, 221.72866928, 222.27632094, 222.8239726,&
223.37162427, 223.91927593, 224.46692759, 225.01457926, 225.56223092,&
226.10988258, 226.65753425, 227.20518591, 227.75283757, 228.30048924,&
228.8481409,  229.39579256, 229.94344423, 230.49109589, 231.03874755,&
231.58639922, 232.13405088, 232.68170254, 233.22935421, 233.77700587,&
234.32465753, 234.8723092,  235.41996086, 235.96761252, 236.51526419,&
237.06291585, 237.61056751, 238.15821918, 238.70587084, 239.2535225,&
239.80117417, 240.34882583, 240.8964775,  241.44412916, 241.99178082,&
242.53943249, 243.08708415, 243.63473581, 244.18238748, 244.73003914,&
245.2776908,  245.82534247, 246.37299413, 246.92064579, 247.46829746,&
248.01594912, 248.56360078, 249.11125245, 249.65890411, 250.20655577,&
250.75420744, 251.3018591,  251.84951076, 252.39716243, 252.94481409,&
253.49246575, 254.04011742, 254.58776908, 255.13542074, 255.68307241,&
256.23072407, 256.77837573, 257.3260274,  257.87367906, 258.42133072,&
258.96898239, 259.51663405, 260.06428571, 260.61193738, 261.15958904,&
261.7072407,  262.25489237, 262.80254403, 263.35019569, 263.89784736,&
264.44549902, 264.99315068, 265.54080235, 266.08845401, 266.63610568,&
267.18375734, 267.731409,   268.27906067, 268.82671233, 269.37436399,&
269.92201566, 270.46966732, 271.01731898, 271.56497065, 272.11262231,&
272.66027397, 273.20792564, 273.7555773,  274.30322896, 274.85088063,&
275.39853229, 275.94618395, 276.49383562, 277.04148728, 277.58913894,&
278.13679061, 278.68444227, 279.23209393, 279.7797456,  280.32739726,&
280.87504892, 281.42270059, 281.97035225, 282.51800391, 283.06565558,&
283.61330724, 284.1609589,  284.70861057, 285.25626223, 285.80391389,&
286.35156556, 286.89921722, 287.44686888, 287.99452055, 288.54217221,&
289.08982387, 289.63747554, 290.1851272,  290.73277886, 291.28043053,&
291.82808219, 292.37573386, 292.92338552, 293.47103718, 294.01868885,&
294.56634051, 295.11399217, 295.66164384, 296.2092955,  296.75694716,&
297.30459883, 297.85225049, 298.39990215, 298.94755382, 299.49520548,&
300.04285714, 300.59050881, 301.13816047, 301.68581213, 302.2334638,&
302.78111546, 303.32876712, 303.87641879, 304.42407045, 304.97172211,&
305.51937378, 306.06702544, 306.6146771,  307.16232877, 307.70998043,&
308.25763209, 308.80528376, 309.35293542, 309.90058708, 310.44823875,&
310.99589041, 311.54354207, 312.09119374, 312.6388454,  313.18649706,&
313.73414873, 314.28180039, 314.82945205, 315.37710372, 315.92475538,&
316.47240705, 317.02005871, 317.56771037, 318.11536204, 318.6630137,&
319.21066536, 319.75831703, 320.30596869, 320.85362035, 321.40127202,&
321.94892368, 322.49657534, 323.04422701, 323.59187867, 324.13953033,&
324.687182  , 325.23483366, 325.78248532, 326.33013699, 326.87778865,&
327.42544031, 327.97309198, 328.52074364, 329.0683953,  329.61604697,&
330.16369863, 330.71135029, 331.25900196, 331.80665362, 332.35430528,&
332.90195695, 333.44960861, 333.99726027, 334.54491194, 335.0925636,&
335.64021526, 336.18786693, 336.73551859, 337.28317025, 337.83082192,&
338.37847358, 338.92612524, 339.47377691, 340.02142857, 340.56908023,&
341.1167319,  341.66438356, 342.21203523, 342.75968689, 343.30733855,&
343.85499022, 344.40264188, 344.95029354, 345.49794521, 346.04559687,&
346.59324853, 347.1409002,  347.68855186, 348.23620352, 348.78385519,&
349.33150685, 349.87915851, 350.42681018, 350.97446184, 351.5221135,&
352.06976517, 352.61741683, 353.16506849, 353.71272016, 354.26037182,&
354.80802348, 355.35567515, 355.90332681, 356.45097847, 356.99863014,&
357.5462818,  358.09393346, 358.64158513, 359.18923679, 359.73688845,&
360.28454012, 360.83219178, 361.37984344, 361.92749511, 362.47514677,&
363.02279843, 363.5704501,  364.11810176, 364.66575342, 365.21340509,&
365.76105675, 366.30870841, 366.85636008, 367.40401174, 367.95166341,&
368.49931507, 369.04696673, 369.5946184,  370.14227006, 370.68992172,&
371.23757339, 371.78522505, 372.33287671, 372.88052838, 373.42818004,&
373.9758317,  374.52348337, 375.07113503, 375.61878669, 376.16643836,&
376.71409002, 377.26174168, 377.80939335, 378.35704501, 378.90469667,&
379.45234834, 380.0 /)

real, parameter :: es_table(512) = (/ 1.00000000e-11,  1.00000000e-11,  1.00000000e-11,  1.00000000e-11,&
   1.00000000e-11,  1.00000000e-11,  1.00000000e-11,  1.00000000e-11,&
   1.00000000e-11, 1.00000000e-11, 1.00000000e-11, 1.00000000e-11,&
   1.00000000e-11, 1.00000000e-11, 1.00000000e-11, 1.00000000e-11,&
   1.00000000e-11, 1.00000000e-11, 1.00000000e-11, 1.00000000e-11,&
   1.00000000e-11, 1.00000000e-11, 1.00000000e-11, 1.04431694e-11,&
   1.35878670e-11, 1.76347924e-11, 2.28299687e-11, 2.94830001e-11,&
   3.79825947e-11, 4.88156526e-11, 6.25906869e-11, 8.00665005e-11,&
   1.02187219e-10, 1.30125003e-10, 1.65332008e-10, 2.09603473e-10,&
   2.65154155e-10, 3.34710762e-10, 4.21623499e-10, 5.30000424e-10,&
   6.64868978e-10, 8.32369796e-10, 1.03998884e-09, 1.29683494e-09,&
   1.61397089e-09, 2.00480804e-09, 2.48557530e-09, 3.07587602e-09,&
   3.79934781e-09, 4.68444322e-09, 5.76535159e-09, 7.08308621e-09,&
   8.68676392e-09, 1.06351092e-08, 1.29982190e-08, 1.58596306e-08,&
   1.93187404e-08, 2.34936295e-08, 2.85243589e-08, 3.45768073e-08,&
   4.18471333e-08, 5.05669578e-08, 6.10093728e-08, 7.34958986e-08,&
   8.84045297e-08, 1.06179025e-07, 1.27339619e-07, 1.52495362e-07,&
   1.82358310e-07, 2.17759819e-07, 2.59669247e-07, 3.09215371e-07,&
   3.67710899e-07, 4.36680483e-07, 5.17892689e-07, 6.13396453e-07,&
   7.25562586e-07, 8.57130984e-07, 1.01126426e-06, 1.19160861e-06,&
   1.40236276e-06, 1.64835609e-06, 1.93513697e-06, 2.26907245e-06,&
   2.65746094e-06, 3.10865906e-06, 3.63222456e-06, 4.23907706e-06,&
   4.94167869e-06, 5.75423689e-06, 6.69293175e-06, 7.77617077e-06,&
   9.02487402e-06, 1.04627928e-05, 1.21168660e-05, 1.40176171e-05,&
   1.61995975e-05, 1.87018798e-05, 2.15686073e-05, 2.48496044e-05,&
   2.86010545e-05, 3.28862535e-05, 3.77764439e-05, 4.33517397e-05,&
   4.97021504e-05, 5.69287125e-05, 6.51447410e-05, 7.44772107e-05,&
   8.50682803e-05, 9.70769720e-05, 1.10681022e-04, 1.26078916e-04,&
   1.43492132e-04, 1.63167594e-04, 1.85380382e-04, 2.10436691e-04,&
   2.38677081e-04, 2.70480039e-04, 3.06265875e-04, 3.46500984e-04,&
   3.91702507e-04, 4.42443421e-04, 4.99358092e-04, 5.63148336e-04,&
   6.34590025e-04, 7.14540280e-04, 8.03945303e-04, 9.03848893e-04,&
   1.01540171e-03, 1.13987133e-03, 1.27865315e-03, 1.43328229e-03,&
   1.60544636e-03, 1.79699944e-03, 2.00997711e-03, 2.24661278e-03,&
   2.50935534e-03, 2.80088817e-03, 3.12414981e-03, 3.48235611e-03,&
   3.87902428e-03, 4.31799870e-03, 4.80347882e-03, 5.34004916e-03,&
   5.93271164e-03, 6.58692029e-03, 7.30861868e-03, 8.10428002e-03,&
   8.98095035e-03, 9.94629478e-03, 1.10086472e-02, 1.21770635e-02,&
   1.34613788e-02, 1.48722684e-02, 1.64213132e-02, 1.81210699e-02,&
   1.99851461e-02, 2.20282800e-02, 2.42664260e-02, 2.67168463e-02,&
   2.93982085e-02, 3.23306891e-02, 3.55360852e-02, 3.90379322e-02,&
   4.28616301e-02, 4.70345775e-02, 5.15863145e-02, 5.65486746e-02,&
   6.19559462e-02, 6.78450445e-02, 7.42556937e-02, 8.12306209e-02,&
   8.88157616e-02, 9.70604783e-02, 1.06017791e-01, 1.15744625e-01,&
   1.26302066e-01, 1.37755642e-01, 1.50175609e-01, 1.63637264e-01,&
   1.78221267e-01, 1.94013993e-01, 2.11107888e-01, 2.29601863e-01,&
   2.49601698e-01, 2.71220468e-01, 2.94579005e-01, 3.19806371e-01,&
   3.47040366e-01, 3.76428056e-01, 4.08126341e-01, 4.42302541e-01,&
   4.79135018e-01, 5.18813836e-01, 5.61541445e-01, 6.07533410e-01,&
   6.57019171e-01, 7.10242848e-01, 7.67464078e-01, 8.28958906e-01,&
   8.95020708e-01, 9.65961172e-01, 1.04211132e+00, 1.12382257e+00,&
   1.21146789e+00, 1.30544294e+00, 1.40616735e+00, 1.51408597e+00,&
   1.62967027e+00, 1.75341975e+00, 1.88586338e+00, 2.02756124e+00,&
   2.17910608e+00, 2.34112503e+00, 2.51428142e+00, 2.69927660e+00,&
   2.89685188e+00, 3.10779057e+00, 3.33292012e+00, 3.57311428e+00,&
   3.82929545e+00, 4.10243705e+00, 4.39356606e+00, 4.70376563e+00,&
   5.03417776e+00, 5.38600621e+00, 5.76051941e+00, 6.15905356e+00,&
   6.58301582e+00, 7.03388768e+00, 7.51322840e+00, 8.02267864e+00,&
   8.56396422e+00, 9.13890000e+00, 9.74939400e+00, 1.03974515e+01,&
   1.10851797e+01, 1.18147918e+01, 1.25886123e+01, 1.34090813e+01,&
   1.42787602e+01, 1.52003366e+01, 1.61766298e+01, 1.72105966e+01,&
   1.83053373e+01, 1.94641013e+01, 2.06902944e+01, 2.19790653e+01,&
   2.33282240e+01, 2.47452937e+01, 2.62345742e+01, 2.77999453e+01,&
   2.94452092e+01, 3.11741989e+01, 3.29908271e+01, 3.48991119e+01,&
   3.69031934e+01, 3.90073458e+01, 4.12159860e+01, 4.35336814e+01,&
   4.59651566e+01, 4.85152999e+01, 5.11891688e+01, 5.39919960e+01,&
   5.69291947e+01, 6.00063645e+01, 6.32292968e+01, 6.66039801e+01,&
   7.01366061e+01, 7.38335748e+01, 7.77015006e+01, 8.17472178e+01,&
   8.59777866e+01, 9.04004990e+01, 9.50228846e+01, 9.98527171e+01,&
   1.04898020e+02, 1.10167073e+02, 1.15668420e+02, 1.21410870e+02,&
   1.27403513e+02, 1.33655717e+02, 1.40177141e+02, 1.46977739e+02,&
   1.54067768e+02, 1.61457795e+02, 1.69158705e+02, 1.77181704e+02,&
   1.85538332e+02, 1.94240468e+02, 2.03300335e+02, 2.12730511e+02,&
   2.22543937e+02, 2.32753921e+02, 2.43374151e+02, 2.54418697e+02,&
   2.65902025e+02, 2.77839001e+02, 2.90244900e+02, 3.03135418e+02,&
   3.16526676e+02, 3.30435227e+02, 3.44878073e+02, 3.59872665e+02,&
   3.75436915e+02, 3.91589208e+02, 4.08348405e+02, 4.25733856e+02,&
   4.43765408e+02, 4.62463417e+02, 4.81848752e+02, 5.01942809e+02,&
   5.22767516e+02, 5.44345350e+02, 5.66699337e+02, 5.89853070e+02,&
   6.13832479e+02, 6.38680126e+02, 6.64428306e+02, 6.91105655e+02,&
   7.18741584e+02, 7.47366296e+02, 7.77010802e+02, 8.07706940e+02,&
   8.39487391e+02, 8.72385701e+02, 9.06436294e+02, 9.41674493e+02,&
   9.78136541e+02, 1.01585962e+03, 1.05488185e+03, 1.09524236e+03,&
   1.13698125e+03, 1.18013965e+03, 1.22475971e+03, 1.27088465e+03,&
   1.31855878e+03, 1.36782749e+03, 1.41873729e+03, 1.47133586e+03,&
   1.52567202e+03, 1.58179580e+03, 1.63975842e+03, 1.69961236e+03,&
   1.76141133e+03, 1.82521036e+03, 1.89106576e+03, 1.95903516e+03,&
   2.02917759e+03, 2.10155342e+03, 2.17622445e+03, 2.25325390e+03,&
   2.33270645e+03, 2.41464828e+03, 2.49914707e+03, 2.58627203e+03,&
   2.67609395e+03, 2.76868522e+03, 2.86411983e+03, 2.96247343e+03,&
   3.06382336e+03, 3.16824865e+03, 3.27583008e+03, 3.38665020e+03,&
   3.50079334e+03, 3.61834566e+03, 3.73939521e+03, 3.86403187e+03,&
   3.99234750e+03, 4.12443587e+03, 4.26039275e+03, 4.40031593e+03,&
   4.54430523e+03, 4.69246258e+03, 4.84489199e+03, 5.00169966e+03,&
   5.16299393e+03, 5.32888538e+03, 5.49948684e+03, 5.67491343e+03,&
   5.85528257e+03, 6.04071407e+03, 6.23133009e+03, 6.42725527e+03,&
   6.62861668e+03, 6.83554391e+03, 7.04816908e+03, 7.26662689e+03,&
   7.49105466e+03, 7.72159237e+03, 7.95838269e+03, 8.20157100e+03,&
   8.45130548e+03, 8.70773712e+03, 8.97101974e+03, 9.24131008e+03,&
   9.51876777e+03, 9.80355546e+03, 1.00958388e+04, 1.03957865e+04,&
   1.07035703e+04, 1.10193651e+04, 1.13433492e+04, 1.16757038e+04,&
   1.20166136e+04, 1.23662664e+04, 1.27248536e+04, 1.30925699e+04,&
   1.34696133e+04, 1.38561856e+04, 1.42524918e+04, 1.46587407e+04,&
   1.50751447e+04, 1.55019198e+04, 1.59392858e+04, 1.63874661e+04,&
   1.68466882e+04, 1.73171832e+04, 1.77991862e+04, 1.82929363e+04,&
   1.87986766e+04, 1.93166542e+04, 1.98471202e+04, 2.03903301e+04,&
   2.09465435e+04, 2.15160241e+04, 2.20990401e+04, 2.26958640e+04,&
   2.33067727e+04, 2.39320474e+04, 2.45719740e+04, 2.52268429e+04,&
   2.58969490e+04, 2.65825922e+04, 2.72840766e+04, 2.80017114e+04,&
   2.87358106e+04, 2.94866930e+04, 3.02546823e+04, 3.10401073e+04,&
   3.18433016e+04, 3.26646041e+04, 3.35043587e+04, 3.43629147e+04,&
   3.52406264e+04, 3.61378536e+04, 3.70549612e+04, 3.79923198e+04,&
   3.89503054e+04, 3.99292994e+04, 4.09296889e+04, 4.19518666e+04,&
   4.29962308e+04, 4.40631858e+04, 4.51531415e+04, 4.62665137e+04,&
   4.74037242e+04, 4.85652007e+04, 4.97513770e+04, 5.09626930e+04,&
   5.21995946e+04, 5.34625342e+04, 5.47519703e+04, 5.60683678e+04,&
   5.74121978e+04, 5.87839381e+04, 6.01840729e+04, 6.16130931e+04,&
   6.30714960e+04, 6.45597857e+04, 6.60784732e+04, 6.76280760e+04,&
   6.92091188e+04, 7.08221331e+04, 7.24676572e+04, 7.41462369e+04,&
   7.58584247e+04, 7.76047804e+04, 7.93858712e+04, 8.12022714e+04,&
   8.30545628e+04, 8.49433345e+04, 8.68691832e+04, 8.88327131e+04,&
   9.08345360e+04, 9.28752715e+04, 9.49555466e+04, 9.70759966e+04,&
   9.92372642e+04, 1.01440000e+05, 1.03684864e+05, 1.05972521e+05,&
   1.08303648e+05, 1.10678927e+05, 1.13099049e+05, 1.15564715e+05,&
   1.18076632e+05, 1.20635517e+05, 1.23242095e+05, 1.25897099e+05,&
   1.28601272e+05, 1.31355364e+05, 1.34160135e+05, 1.37016354e+05,&
   1.39924797e+05, 1.42886251e+05, 1.45901511e+05, 1.48971382e+05,&
   1.52096676e+05, 1.55278216e+05, 1.58516835e+05, 1.61813373e+05/)

  real, parameter :: des_table(512) = (/6.12116769e-12, 6.05476789e-12, 5.98944267e-12, 5.92516898e-12,&
   5.86192436e-12, 5.79968696e-12, 5.73843551e-12, 5.67814929e-12,&
   5.61880812e-12, 5.56039236e-12, 5.50288287e-12, 5.44626098e-12,&
   5.39050854e-12, 5.33560783e-12, 5.28154159e-12, 5.22829299e-12,&
   5.17584563e-12, 5.12418352e-12, 5.07329105e-12, 5.02315302e-12,&
   4.97375458e-12, 4.92508127e-12, 4.87711895e-12, 5.04389820e-12,&
   6.49944708e-12, 8.35423924e-12, 1.07120700e-11, 1.37022285e-11,&
   1.74853900e-11, 2.22607857e-11, 2.82749050e-11, 3.58320369e-11,&
   4.53070134e-11, 5.71605796e-11, 7.19578978e-11, 9.03907779e-11,&
   1.13304331e-10, 1.41728863e-10, 1.76917973e-10, 2.20393961e-10,&
   2.74001862e-10, 3.39973616e-10, 4.21004144e-10, 5.20341371e-10,&
   6.41892564e-10, 7.90349728e-10, 9.71337222e-10, 1.19158522e-09,&
   1.45913324e-09, 1.78356846e-09, 2.17630454e-09, 2.65090694e-09,&
   3.22347225e-09, 3.91306966e-09, 4.74225383e-09, 5.73766007e-09,&
   6.93069370e-09, 8.35832761e-09, 1.00640234e-08, 1.20987938e-08,&
   1.45224265e-08, 1.74048914e-08, 2.08279568e-08, 2.48870441e-08,&
   2.96933505e-08, 3.53762784e-08, 4.20862101e-08, 4.99976735e-08,&
   5.93129489e-08, 7.02661735e-08, 8.31280066e-08, 9.82109238e-08,&
   1.15875221e-07, 1.36535812e-07, 1.60669915e-07, 1.88825741e-07,&
   2.21632290e-07, 2.59810396e-07, 3.04185164e-07, 3.55699948e-07,&
   4.15432058e-07, 4.84610382e-07, 5.64635133e-07, 6.57099966e-07,&
   7.63816706e-07, 8.86842982e-07, 1.02851306e-06, 1.19147224e-06,&
   1.37871514e-06, 1.59362832e-06, 1.84003766e-06, 2.12226095e-06,&
   2.44516629e-06, 2.81423677e-06, 3.23564212e-06, 3.71631797e-06,&
   4.26405349e-06, 4.88758811e-06, 5.59671830e-06, 6.40241526e-06,&
   7.31695450e-06, 8.35405851e-06, 9.52905357e-06, 1.08590420e-05,&
   1.23630913e-05, 1.40624415e-05, 1.59807323e-05, 1.81442524e-05,&
   2.05822110e-05, 2.33270359e-05, 2.64146983e-05, 2.98850680e-05,&
   3.37822998e-05, 3.81552556e-05, 4.30579636e-05, 4.85501178e-05,&
   5.46976216e-05, 6.15731780e-05, 6.92569309e-05, 7.78371604e-05,&
   8.74110374e-05, 9.80854399e-05, 1.09977839e-04, 1.23217254e-04,&
   1.37945290e-04, 1.54317255e-04, 1.72503368e-04, 1.92690062e-04,&
   2.15081387e-04, 2.39900530e-04, 2.67391443e-04, 2.97820598e-04,&
   3.31478876e-04, 3.68683602e-04, 4.09780718e-04, 4.55147134e-04,&
   5.05193234e-04, 5.60365577e-04, 6.21149781e-04, 6.88073623e-04,&
   7.61710344e-04, 8.42682203e-04, 9.31664258e-04, 1.02938842e-03,&
   1.13664779e-03, 1.25430126e-03, 1.38327846e-03, 1.52458499e-03,&
   1.67930804e-03, 1.84862237e-03, 2.03379663e-03, 2.23620016e-03,&
   2.45731013e-03, 2.69871925e-03, 2.96214385e-03, 3.24943251e-03,&
   3.56257523e-03, 3.90371316e-03, 4.27514887e-03, 4.67935730e-03,&
   5.11899733e-03, 5.59692399e-03, 6.11620152e-03, 6.68011698e-03,&
   7.29219487e-03, 7.95621238e-03, 8.67621565e-03, 9.45653685e-03,&
   1.03018123e-02, 1.12170014e-02, 1.22074069e-02, 1.32786960e-02,&
   1.44369227e-02, 1.56885513e-02, 1.70404812e-02, 1.85000729e-02,&
   2.00751756e-02, 2.17741557e-02, 2.36059277e-02, 2.55799852e-02,&
   2.77064351e-02, 2.99960325e-02, 3.24602173e-02, 3.51111532e-02,&
   3.79617681e-02, 4.10257968e-02, 4.43178250e-02, 4.78533369e-02,&
   5.16487630e-02, 5.57215321e-02, 6.00901243e-02, 6.47741274e-02,&
   6.97942949e-02, 7.51726076e-02, 8.09323376e-02, 8.70981146e-02,&
   9.36959959e-02, 1.00753539e-01, 1.08299879e-01, 1.16365804e-01,&
   1.24983843e-01, 1.34188347e-01, 1.44015581e-01, 1.54503818e-01,&
   1.65693437e-01, 1.77627020e-01, 1.90349463e-01, 2.03908080e-01,&
   2.18352723e-01, 2.33735897e-01, 2.50112882e-01, 2.67541865e-01,&
   2.86084071e-01, 3.05803901e-01, 3.26769075e-01, 3.49050782e-01,&
   3.72723835e-01, 3.97866828e-01, 4.24562307e-01, 4.52896937e-01,&
   4.82961685e-01, 5.14852002e-01, 5.48668016e-01, 5.84514728e-01,&
   6.22502221e-01, 6.62745869e-01, 7.05366559e-01, 7.50490917e-01,&
   7.98251544e-01, 8.48787257e-01, 9.02243344e-01, 9.58771820e-01,&
   1.01853170e+00, 1.08168926e+00, 1.14841835e+00, 1.21890066e+00,&
   1.29332605e+00, 1.37189285e+00, 1.45480819e+00, 1.54228831e+00,&
   1.63455894e+00, 1.73185565e+00, 1.83442418e+00, 1.94252086e+00,&
   2.05641297e+00, 2.17637916e+00, 2.30270985e+00, 2.40474715e+00,&
   2.52409614e+00, 2.65225593e+00, 2.78770444e+00, 2.93011692e+00,&
   3.07948018e+00, 3.23590497e+00, 3.39956432e+00, 3.57066807e+00,&
   3.74945089e+00, 3.93616595e+00, 4.13108146e+00, 4.33447863e+00,&
   4.54665040e+00, 4.76790075e+00, 4.99854421e+00, 5.23890562e+00,&
   5.48931998e+00, 5.75013240e+00, 6.02169810e+00, 6.30438248e+00,&
   6.59856117e+00, 6.90462016e+00, 7.22295588e+00, 7.55397538e+00,&
   7.89809642e+00, 8.25574766e+00, 8.62736877e+00, 9.01341065e+00,&
   9.41433552e+00, 9.83061713e+00, 1.02627409e+01, 1.07112042e+01,&
   1.11765164e+01, 1.16591988e+01, 1.21597856e+01, 1.26788231e+01,&
   1.32168705e+01, 1.37744999e+01, 1.43522965e+01, 1.49508586e+01,&
   1.55707983e+01, 1.62127408e+01, 1.68773255e+01, 1.75652055e+01,&
   1.82770481e+01, 1.90135350e+01, 1.97753623e+01, 2.05632406e+01,&
   2.13778957e+01, 2.22200680e+01, 2.30905133e+01, 2.39900026e+01,&
   2.49193225e+01, 2.58792753e+01, 2.68706789e+01, 2.78943676e+01,&
   2.89511914e+01, 3.00420169e+01, 3.11677271e+01, 3.23292218e+01,&
   3.35274174e+01, 3.47632472e+01, 3.60376619e+01, 3.73516292e+01,&
   3.87061344e+01, 4.01021800e+01, 4.15407867e+01, 4.30229926e+01,&
   4.45661649e+01, 4.61848393e+01, 4.78551049e+01, 4.95783612e+01,&
   5.13560386e+01, 5.31895982e+01, 5.50805333e+01, 5.70303689e+01,&
   5.90406628e+01, 6.11130061e+01, 6.32490233e+01, 6.54503732e+01,&
   6.77187492e+01, 7.00558799e+01, 7.24635298e+01, 7.49434992e+01,&
   7.74976256e+01, 8.01277835e+01, 8.28358854e+01, 8.56238819e+01,&
   8.84937628e+01, 9.14475572e+01, 9.44873340e+01, 9.76152030e+01,&
   1.00833315e+02, 1.04143862e+02, 1.07549078e+02, 1.11051241e+02,&
   1.14652672e+02, 1.18355734e+02, 1.22162837e+02, 1.26076434e+02,&
   1.30099026e+02, 1.34233157e+02, 1.38481420e+02, 1.42846455e+02,&
   1.47330949e+02, 1.51937639e+02, 1.56669309e+02, 1.61528794e+02,&
   1.66518980e+02, 1.71642801e+02, 1.76903247e+02, 1.82303355e+02,&
   1.87846218e+02, 1.93534981e+02, 1.99372842e+02, 2.05363057e+02,&
   2.11508932e+02, 2.17813832e+02, 2.24281177e+02, 2.30914444e+02,&
   2.37717169e+02, 2.44692943e+02, 2.51845418e+02, 2.59178305e+02,&
   2.66695373e+02, 2.74400455e+02, 2.82297441e+02, 2.90390287e+02,&
   2.98683007e+02, 3.07179681e+02, 3.15884452e+02, 3.24801527e+02,&
   3.33935176e+02, 3.43289737e+02, 3.52869613e+02, 3.62679274e+02,&
   3.72723256e+02, 3.83006164e+02, 3.93532671e+02, 4.04307520e+02,&
   4.15335522e+02, 4.26621560e+02, 4.38170586e+02, 4.49987625e+02,&
   4.62077774e+02, 4.74446203e+02, 4.87098153e+02, 5.00038943e+02,&
   5.13273963e+02, 5.26808679e+02, 5.40648633e+02, 5.54799445e+02,&
   5.69266810e+02, 5.84056500e+02, 5.99174366e+02, 6.14626338e+02,&
   6.30418427e+02, 6.46556719e+02, 6.63047385e+02, 6.79896676e+02,&
   6.97110923e+02, 7.14696542e+02, 7.32660028e+02, 7.51007964e+02,&
   7.69747013e+02, 7.88883924e+02, 8.08425532e+02, 8.28378757e+02,&
   8.48750604e+02, 8.69548166e+02, 8.90778624e+02, 9.12449245e+02,&
   9.34567386e+02, 9.57140493e+02, 9.80176099e+02, 1.00368183e+03,&
   1.02766540e+03, 1.05213462e+03, 1.07709739e+03, 1.10256169e+03,&
   1.12853561e+03, 1.15502733e+03, 1.18204511e+03, 1.20959732e+03,&
   1.23769243e+03, 1.26633897e+03, 1.29554561e+03, 1.32532109e+03,&
   1.35567425e+03, 1.38661403e+03, 1.41814947e+03, 1.45028971e+03,&
   1.48304398e+03, 1.51642162e+03, 1.55043205e+03, 1.58508482e+03,&
   1.62038954e+03, 1.65635597e+03, 1.69299393e+03, 1.73031336e+03,&
   1.76832431e+03, 1.80703690e+03, 1.84646140e+03, 1.88660814e+03,&
   1.92748759e+03, 1.96911029e+03, 2.01148691e+03, 2.05462822e+03,&
   2.09854508e+03, 2.14324848e+03, 2.18874949e+03, 2.23505932e+03,&
   2.28218925e+03, 2.33015070e+03, 2.37895516e+03, 2.42861427e+03,&
   2.47913976e+03, 2.53054345e+03, 2.58283730e+03, 2.63603336e+03,&
   2.69014380e+03, 2.74518089e+03, 2.80115702e+03, 2.85808467e+03,&
   2.91597647e+03, 2.97484512e+03, 3.03470346e+03, 3.09556443e+03,&
   3.15744107e+03, 3.22034655e+03, 3.28429415e+03, 3.34929727e+03,&
   3.41536939e+03, 3.48252414e+03, 3.55077525e+03, 3.62013656e+03,&
   3.69062202e+03, 3.76224571e+03, 3.83502181e+03, 3.90896461e+03,&
   3.98408854e+03, 4.06040811e+03, 4.13793798e+03, 4.21669290e+03,&
   4.29668774e+03, 4.37793750e+03, 4.46045726e+03, 4.54426227e+03,&
   4.62936784e+03, 4.71578944e+03, 4.80354262e+03, 4.89264307e+03,&
   4.98310660e+03, 5.07494911e+03, 5.16818664e+03, 5.26283533e+03,&
   5.35891146e+03, 5.45643140e+03, 5.55541165e+03, 5.65586883e+03,&
   5.75781966e+03, 5.86128100e+03, 5.96626981e+03, 6.07280317e+03/)

  real, parameter :: dx = t_table(2) - t_table(1)
  real, parameter :: dxi = 1.0 / dx
  real, parameter :: x_min = t_table(1)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


contains

!#######################################################################

 subroutine lookup_es_0d ( temp, esat )

 real, intent(in)  :: temp
 real, intent(out) :: esat
 integer :: i, ind
 real :: y1, y2, x1
!-----------------------------------------------

! ind = minloc(abs(t_table - temp), 1)
! esat = es_table(ind)
ind = floor((temp - x_min)*dxi)
y1 = es_table(ind)
y2 = es_table(ind + 1)
x1 = x_min + dx * ind
esat = y1 + (temp - x1) * (y2 - y1) * dxi

!-----------------------------------------------

 end subroutine lookup_es_0d

!#######################################################################

 subroutine lookup_es_1d ( temp, esat )

 real, intent(in)  :: temp(:)
 real, intent(out) :: esat(:)
 integer :: i, ind
 real :: y1, y2, x1
!-----------------------------------------------

    do i = 1, size(temp)
      ! ind = minloc(abs(t_table - temp(i)), 1)
      ! esat(i) = es_table(ind)
      ind = floor((temp(i) - x_min)*dxi)
      y1 = es_table(ind)
      y2 = es_table(ind + 1)
      x1 = x_min + dx * ind
      esat(i) = y1 + (temp(i) - x1) * (y2 - y1) * dxi
    enddo
!-----------------------------------------------

 end subroutine lookup_es_1d

!#######################################################################

 subroutine lookup_es_2d ( temp, esat )

 real, intent(in)  :: temp(:,:)
 real, intent(out) :: esat(:,:)
 integer :: i, j, ind
 real :: y1, y2, x1
!-----------------------------------------------

    do i = 1, size(temp,1)
    do j = 1, size(temp,2)
      ! ind = minloc(abs(t_table - temp(i, j)), 1)
      ! esat(i,j) = es_table(ind)
      ind = floor((temp(i, j) - x_min)*dxi)
      y1 = es_table(ind)
      y2 = es_table(ind + 1)
      x1 = x_min + dx * ind
      esat(i, j) = y1 + (temp(i, j) - x1) * (y2 - y1) * dxi
    enddo
    enddo

!-----------------------------------------------

 end subroutine lookup_es_2d

!#######################################################################

 subroutine lookup_es_3d ( temp, esat )

 real, intent(in)  :: temp(:,:,:)
 real, intent(out) :: esat(:,:,:)
 integer :: i, j, k, ind
 real :: y1, y2, x1
!-----------------------------------------------

    do i = 1, size(temp,1)
    do j = 1, size(temp,2)
    do k = 1, size(temp,3)
      ! ind = minloc(abs(t_table - temp(i,j,k)), 1)
      ! esat(i,j,k) = es_table(ind)
      ind = floor((temp(i, j, k) - x_min)*dxi)
      y1 = es_table(ind)
      y2 = es_table(ind + 1)
      x1 = x_min + dx * ind
      esat(i, j, k) = y1 + (temp(i, j, k) - x1) * (y2 - y1) * dxi
    enddo
    enddo
    enddo
!-----------------------------------------------

 end subroutine lookup_es_3d

 !#######################################################################

  subroutine lookup_des_0d ( temp, desat )

  real, intent(in)  :: temp
  real, intent(out) :: desat
  integer :: i, ind
  real :: y1, y2, x1
 !-----------------------------------------------

 ind = floor((temp - x_min)*dxi)
 y1 = des_table(ind)
 y2 = des_table(ind + 1)
 x1 = x_min + dx * ind
 desat = y1 + (temp - x1) * (y2 - y1) * dxi

 !-----------------------------------------------

  end subroutine lookup_des_0d

 !#######################################################################

  subroutine lookup_des_1d ( temp, desat )

  real, intent(in)  :: temp(:)
  real, intent(out) :: desat(:)
  integer :: i, ind
  real :: y1, y2, x1
 !-----------------------------------------------

     do i = 1, size(temp)
       ind = floor((temp(i) - x_min)*dxi)
       y1 = des_table(ind)
       y2 = des_table(ind + 1)
       x1 = x_min + dx * ind
       desat(i) = y1 + (temp(i) - x1) * (y2 - y1) * dxi
     enddo
 !-----------------------------------------------

  end subroutine lookup_des_1d

 !#######################################################################

  subroutine lookup_des_2d ( temp, desat )

  real, intent(in)  :: temp(:,:)
  real, intent(out) :: desat(:,:)
  integer :: i, j, ind
  real :: y1, y2, x1
 !-----------------------------------------------

     do i = 1, size(temp,1)
     do j = 1, size(temp,2)
       ind = floor((temp(i, j) - x_min)*dxi)
       y1 = des_table(ind)
       y2 = des_table(ind + 1)
       x1 = x_min + dx * ind
       desat(i, j) = y1 + (temp(i, j) - x1) * (y2 - y1) * dxi
     enddo
     enddo

 !-----------------------------------------------

  end subroutine lookup_des_2d

 !#######################################################################

  subroutine lookup_des_3d ( temp, desat )

  real, intent(in)  :: temp(:,:,:)
  real, intent(out) :: desat(:,:,:)
  integer :: i, j, k, ind
  real :: y1, y2, x1
 !-----------------------------------------------

     do i = 1, size(temp,1)
     do j = 1, size(temp,2)
     do k = 1, size(temp,3)
       ind = floor((temp(i, j, k) - x_min)*dxi)
       y1 = des_table(ind)
       y2 = des_table(ind + 1)
       x1 = x_min + dx * ind
       desat(i, j, k) = y1 + (temp(i, j, k) - x1) * (y2 - y1) * dxi
     enddo
     enddo
     enddo
 !-----------------------------------------------

  end subroutine lookup_des_3d


end module simple_sat_vapor_pres_mod
