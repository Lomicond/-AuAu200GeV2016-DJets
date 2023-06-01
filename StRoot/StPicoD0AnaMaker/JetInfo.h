/*
Some additional info and functions for the analysis
*/


#include "StPicoEvent/StPicoDst.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoBEmcPidTraits.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StPicoCuts/StPicoCuts.h"

#include "StEmcUtil/geometry/StEmcGeom.h"


//#include "../StPicoJetMaker/StPicoJetMaker.h"
#include "../StRefMultCorr/StRefMultCorr.h"

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TPythia6.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include "TDatime.h"

#include <../../fastjet/config.h>
#include <../../fastjet/PseudoJet.hh>
#include <../../fastjet/JetDefinition.hh>
#include <../../fastjet/ClusterSequence.hh>
#include <../../fastjet/ClusterSequenceArea.hh>


#include <vector>

using namespace std;
using namespace fastjet;

const int BadTowerArrT[433]={31, 34, 38, 59, 95, 106, 113, 134, 139, 157, 
193, 214, 257, 266, 267, 282, 286, 287, 293, 371, 
385, 405, 410, 426, 433, 460, 474, 504, 506, 533, 
541, 555, 560, 561, 562, 615, 616, 633, 637, 638, 
650, 653, 657, 673, 693, 740, 749, 757, 758, 779, 
789, 790, 791, 792, 793, 796, 799, 803, 806, 809, 
810, 811, 812, 813, 814, 817, 821, 822, 823, 824, 
829, 830, 831, 832, 835, 837, 841, 842, 843, 844, 
846, 849, 850, 851, 852, 853, 857, 873, 875, 893, 
897, 899, 903, 916, 924, 939, 953, 954, 956, 989, 
993, 1005, 1012, 1020, 1023, 1026, 1027, 1028, 1039, 1040, 
1042, 1044, 1045, 1046, 1048, 1057, 1080, 1081, 1100, 1125, 
1130, 1132, 1154, 1159, 1160, 1165, 1171, 1180, 1187, 1189, 
1190, 1197, 1198, 1199, 1200, 1202, 1207, 1208, 1214, 1217, 
1218, 1219, 1220, 1221, 1222, 1223, 1224, 1237, 1238, 1240, 
1241, 1242, 1243, 1244, 1257, 1258, 1259, 1260, 1274, 1284, 
1293, 1298, 1304, 1312, 1329, 1337, 1341, 1348, 1353, 1354, 
1369, 1375, 1378, 1382, 1388, 1394, 1401, 1407, 1408, 1409, 
1427, 1434, 1440, 1448, 1475, 1486, 1487, 1537, 1567, 1574, 
1575, 1588, 1597, 1599, 1612, 1654, 1668, 1679, 1701, 1702, 
1705, 1709, 1713, 1720, 1728, 1745, 1762, 1765, 1766, 1781, 
1786, 1789, 1807, 1819, 1823, 1856, 1866, 1877, 1878, 1901, 
1945, 1984, 2032, 2040, 2073, 2077, 2092, 2093, 2097, 2104, 
2128, 2129, 2168, 2192, 2196, 2214, 2222, 2223, 2278, 2290, 
2303, 2309, 2310, 2311, 2312, 2366, 2386, 2390, 2391, 2392, 
2409, 2415, 2417, 2439, 2458, 2459, 2497, 2521, 2589, 2590, 
2711, 2749, 2781, 2782, 2783, 2784, 2801, 2802, 2803, 2804, 
2816, 2821, 2822, 2823, 2824, 2834, 2841, 2842, 2843, 2844, 
2865, 2890, 2929, 2961, 2969, 2973, 2974, 2975, 2976, 2977,
2978, 2994, 3005, 3017, 3028, 3045, 3056, 3070, 3071, 3146,
3186, 3220, 3263, 3299, 3316, 3320, 3328, 3329, 3337, 3349, 
3350, 3351, 3352, 3354, 3355, 3356, 3360, 3362, 3369, 3370, 
3371, 3372, 3385, 3386, 3397, 3405, 3407, 3417, 3418, 3419, 
3425, 3426, 3427, 3428, 3432, 3433, 3434, 3435, 3436, 3438, 
3445, 3446, 3447, 3448, 3452, 3454, 3455, 3456, 3469, 3487, 
3493, 3494, 3495, 3498, 3499, 3508, 3514, 3516, 3584, 3588, 
3594, 3595, 3599, 3604, 3611, 3668, 3670, 3678, 3679, 3690, 
3692, 3717, 3718, 3720, 3725, 3726, 3732, 3738, 3739, 3757, 
3769, 3777, 3822, 3838, 3840, 3984, 4006, 4013, 4017, 4018, 
4019, 4053, 4057, 4059, 4099, 4139, 4175, 4195, 4198, 4199, 
4217, 4218, 4223, 4238, 4299, 4312, 4331, 4339, 4350, 4355, 
4357, 4405, 4438, 4457, 4458, 4469, 4495, 4496, 4497, 4498, 
4499, 4500, 4539, 4558, 4560, 4653, 4677, 4678, 4684, 4763, 
4768, 4778, 4783}; //Hanseul's bad run list - 403 - with tight dead tower selection 7/31/19 + 32 Dan's swapped towers - 2 duplicate

const int BadTowerArr[822]={31, 34, 35, 38, 95, 96, 106, 113, 114, 134, 139, 157, 160, 175, 193, 200, 214, 220, 224, 257, 
266, 267, 282, 286, 287, 296, 315, 317, 319, 340, 365, 371, 380, 389, 395, 405, 410, 420, 426, 433, 
474, 483, 484, 504, 506, 520, 529, 533, 541, 555, 560, 561, 562, 580, 582, 584, 585, 600, 615, 617, 
633, 635, 637, 638, 643, 649, 650, 653, 657, 671, 673, 674, 677, 693, 708, 749, 757, 759, 775, 776, 
779, 783, 790, 793, 796, 799, 803, 810, 812, 813, 814, 817, 822, 825, 832, 837, 840, 844, 846, 848, 
853, 857, 859, 873, 875, 887, 893, 897, 899, 903, 916, 924, 939, 940, 946, 953, 954, 956, 972, 979, 
980, 989, 993, 996, 997, 999, 1005, 1012, 1014, 1016, 1017, 1020, 1023, 1026, 1027, 1028, 1039, 1040, 1042, 1044, 
1045, 1046, 1048, 1055, 1056, 1057, 1059, 1062, 1064, 1078, 1080, 1081, 1083, 1084, 1090, 1100, 1101, 1102, 1103, 1104, 
1122, 1124, 1125, 1127, 1128, 1130, 1132, 1137, 1141, 1142, 1143, 1144, 1145, 1146, 1147, 1148, 1149, 1150, 1151, 1152, 
1153, 1154, 1155, 1156, 1157, 1158, 1159, 1160, 1161, 1162, 1163, 1164, 1165, 1166, 1167, 1168, 1169, 1170, 1171, 1172, 
1173, 1174, 1175, 1176, 1177, 1178, 1179, 1180, 1181, 1182, 1183, 1184, 1185, 1186, 1187, 1188, 1189, 1190, 1191, 1192, 
1193, 1194, 1195, 1196, 1197, 1198, 1199, 1200, 1201, 1202, 1203, 1204, 1205, 1206, 1207, 1208, 1209, 1210, 1211, 1212, 
1213, 1214, 1215, 1216, 1217, 1218, 1219, 1220, 1221, 1224, 1232, 1237, 1238, 1239, 1240, 1244, 1250, 1257, 1258, 1259, 
1260, 1262, 1274, 1279, 1280, 1284, 1288, 1293, 1294, 1298, 1300, 1304, 1307, 1308, 1312, 1313, 1325, 1329, 1335, 1337, 
1340, 1341, 1348, 1353, 1354, 1366, 1369, 1375, 1376, 1378, 1388, 1394, 1405, 1407, 1408, 1409, 1434, 1436, 1439, 1440, 
1448, 1475, 1480, 1486, 1537, 1567, 1574, 1588, 1592, 1597, 1599, 1612, 1619, 1620, 1654, 1668, 1679, 1701, 1702, 1705, 
1720, 1728, 1740, 1745, 1753, 1759, 1762, 1765, 1766, 1773, 1781, 1786, 1789, 1807, 1856, 1860, 1866, 1877, 1878, 1879, 
1901, 1920, 1938, 1945, 1984, 2000, 2032, 2040, 2059, 2073, 2077, 2080, 2092, 2093, 2097, 2104, 2120, 2128, 2129, 2140, 
2160, 2162, 2168, 2175, 2176, 2177, 2192, 2195, 2196, 2197, 2200, 2202, 2214, 2215, 2216, 2217, 2222, 2223, 2240, 2243, 
2260, 2278, 2299, 2303, 2305, 2309, 2310, 2311, 2312, 2339, 2340, 2357, 2366, 2386, 2390, 2391, 2392, 2409, 2415, 2417, 
2419, 2420, 2439, 2445, 2458, 2459, 2478, 2479, 2497, 2500, 2535, 2539, 2540, 2559, 2560, 2579, 2580, 2582, 2589, 2590, 
2591, 2592, 2596, 2598, 2609, 2610, 2611, 2612, 2619, 2629, 2630, 2631, 2632, 2637, 2639, 2649, 2650, 2651, 2652, 2669, 
2670, 2671, 2672, 2678, 2689, 2690, 2691, 2692, 2709, 2710, 2711, 2712, 2715, 2717, 2718, 2719, 2729, 2730, 2731, 2732, 
2738, 2749, 2753, 2754, 2755, 2756, 2773, 2774, 2775, 2776, 2781, 2782, 2782, 2783, 2784, 2793, 2794, 2795, 2796, 2801, 
2802, 2803, 2804, 2813, 2814, 2815, 2816, 2820, 2821, 2822, 2822, 2823, 2824, 2834, 2835, 2836, 2841, 2842, 2843, 2844, 
2858, 2865, 2874, 2880, 2890, 2918, 2929, 2961, 2969, 2973, 2974, 2975, 2976, 2977, 2978, 2981, 2982, 2983, 2984, 2985, 
2986, 2987, 2988, 2989, 2990, 2991, 2992, 2993, 2994, 2995, 2996, 2997, 2998, 2999, 3000, 3001, 3002, 3003, 3004, 3005, 
3006, 3007, 3008, 3009, 3010, 3011, 3012, 3013, 3014, 3015, 3016, 3017, 3018, 3019, 3020, 3021, 3022, 3023, 3024, 3025, 
3026, 3027, 3028, 3029, 3030, 3031, 3032, 3033, 3034, 3035, 3036, 3037, 3038, 3039, 3040, 3041, 3042, 3043, 3044, 3045, 
3046, 3047, 3048, 3049, 3050, 3051, 3052, 3053, 3054, 3055, 3056, 3057, 3058, 3059, 3060, 3070, 3071, 3079, 3098, 3099, 
3100, 3139, 3146, 3186, 3218, 3220, 3240, 3263, 3288, 3298, 3299, 3300, 3316, 3320, 3328, 3329, 3337, 3339, 3349, 3350, 
3351, 3352, 3354, 3355, 3356, 3360, 3362, 3369, 3370, 3371, 3372, 3377, 3378, 3379, 3380, 3381, 3382, 3383, 3384, 3385, 
3386, 3387, 3388, 3397, 3399, 3403, 3405, 3410, 3417, 3418, 3419, 3420, 3425, 3425, 3426, 3426, 3427, 3427, 3428, 3428, 
3432, 3433, 3433, 3434, 3434, 3435, 3435, 3436, 3436, 3438, 3445, 3445, 3446, 3446, 3447, 3447, 3448, 3448, 3452, 3452, 
3454, 3454, 3455, 3455, 3456, 3456, 3469, 3473, 3479, 3487, 3493, 3494, 3495, 3498, 3499, 3514, 3516, 3518, 3534, 3555, 
3580, 3584, 3588, 3589, 3594, 3595, 3596, 3599, 3600, 3603, 3611, 3616, 3668, 3670, 3678, 3679, 3690, 3692, 3700, 3717, 
3718, 3720, 3725, 3738, 3739, 3757, 3769, 3777, 3780, 3800, 3838, 3840, 3880, 3897, 3984, 4006, 4013, 4017, 4018, 4019, 
4020, 4037, 4038, 4039, 4040, 4053, 4057, 4058, 4059, 4060, 4077, 4078, 4079, 4080, 4097, 4098, 4099, 4100, 4117, 4118, 
4119, 4120, 4124, 4137, 4138, 4139, 4140, 4157, 4158, 4159, 4160, 4171, 4175, 4177, 4178, 4179, 4180, 4217, 4220, 4223, 
4259, 4279, 4288, 4300, 4312, 4318, 4331, 4350, 4355, 4357, 4369, 4400, 4405, 4437, 4438, 4458, 4459, 4469, 4479, 4495, 
4496, 4497, 4498, 4499, 4500, 4519, 4520, 4539, 4557, 4560, 4563, 4579, 4618, 4657, 4677, 4678, 4684, 4717, 4737, 4763, 
4768, 4783}; //Hanseul's bad tower list https://drupal.star.bnl.gov/STAR/blog/hanseul/run-14-bemc-tower-selection-new-calibration


 const int BadTowerMap[4800]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,
 0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,
 0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,
 0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,1,0,0,0,0,1,0,0,0,0,0,1,1,
 0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,1,
 0,0,0,0,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,1,1,1,0,0,1,0,0,0,0,1,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,0,1,0,0,0,1,0,1,0,1,0,0,0,0,1,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,1,0,1,0,0,0,1,0,0,0,0,
 0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,1,1,0,1,0,0,0,0,0,1,0,0,0,0,0,0,1,0,1,0,1,1,0,0,1,0,0,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,1,
 1,0,1,0,1,1,1,0,1,0,0,0,0,0,0,1,1,1,0,1,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,1,1,0,1,0,1,0,0,0,0,1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,1,1,1,1,0,0,0,1,0,
 0,0,0,0,1,0,0,0,0,0,0,1,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,1,0,0,0,1,0,0,0,0,1,1,0,0,0,1,0,1,0,0,0,1,0,0,1,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,1,0,1,0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,
 1,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,1,0,0,0,0,
 0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,
 0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,1,0,0,0,0,
 0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,1,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,
 0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
 0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
 0,1,1,0,1,0,0,0,0,0,0,1,1,1,1,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,
 0,1,0,1,1,1,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,0,
 0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,1,1,1,1,1,1,0,0,1,1,1,1,1,
 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,1,1,1,1,0,1,1,1,0,0,0,1,0,1,0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,
 0,0,1,0,1,0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,1,1,1,1,1,0,1,0,0,0,0,0,0,1,1,1,1,0,0,0,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,1,1,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
 0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,1,0,0,0,0,1,1,1,0,0,1,1,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,
 0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,
 1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,
 0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,
 0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
 0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,
 0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; //based on BadTowerArr[822]



vector<float> TowArr;
vector<float> Ecorr;
vector<float> TowEta;
vector<float> TowPhi;

//tower based
//array<int, 4800> Nmatch;
//array<float, 4800> SumpT;
array<float, 4800> Sump;
//array<float, 4800> E1;
//array<float, 4800> E2;

//cluster based
//array<int, 1200> Nmatch;
//array<float, 1200> SumpT;
//array<float, 1200> Sump;
//array<float, 1200> EneCl;

vector<vector<int> > Clusters; //dont forget to clear it in the source code
array<int, 4800> Nmatched = {0};
array<double, 4800> SumpCl = {0};
vector<float> EnergyCl;

int nBadTowers = sizeof(BadTowerArr)/sizeof(BadTowerArr[0]);



int indexarr[1200][4] = {0};

TPythia6 *fpythia;


/*class JetInfo : public PseudoJet::UserInfoBase {
public:
		bool isbemctrack;
    short index;
		StPicoBEmcPidTraits* emctraits;
		int towid;
		float trkE;
		float pTBEMC;

    //AddJetInfo(float vx, float vy, float vz): vertex(vx, vy, vz) {
    //}
		JetInfo(bool _isbemctrack, short _index, StPicoBEmcPidTraits* _emctraits, int _towid, float _trkE, float _pTBEMC): isbemctrack(_isbemctrack), index(_index), emctraits(_emctraits), towid(_towid), trkE(_trkE), pTBEMC(_pTBEMC) {}
};
*/

double Median(vector<PseudoJet> jets_used) {
						
						double median = 0;
						vector<double> p_T; 	

						if (jets_used.size() > 0) {
						for (unsigned int i = 0; i < jets_used.size(); i++)
						{
						p_T.push_back(jets_used[i].perp());
						p_T[i] /= jets_used[i].area();
						jets_used[i].reset_momentum(p_T[i]/sqrt(2), p_T[i]/sqrt(2), jets_used[i].pz(), jets_used[i].E());
						} 
						p_T.clear();
						vector<PseudoJet> sorted = sorted_by_pt(jets_used);
						if (sorted.size() % 2 != 0)
						{median = sorted[ceil((double)sorted.size()/2)].perp();}
  					else {median = (sorted[sorted.size()/2].perp()+sorted[sorted.size()/2+1].perp())/2;}
						}
	
  return median;
}

//----------------------------------------------------------------------------- 
//Generate jet for embedding
//jetType: 0...single particle jet | 1...Pythia jet (TBD)
//----------------------------------------------------------------------------- 
vector<PseudoJet> EmbedJet(short jetType, float &pT_emb, float maxRapJet /*was float R and maxRapJet was calculated inside*/, vector<PseudoJet> container/*, float* eta_emb, float* phi_emb*/, TString parton, TPythia6* pythia)
{
	//float maxRapJet=fMaxRap - R; //fiducial jet acceptance
	//float maxRapJet=fMaxRap - 0.3; //fiducial jet acceptance
	double eta_rnd = gRandom->Uniform(-(maxRapJet),(maxRapJet));
   double phi_rnd = gRandom->Uniform(0, 2.*TMath::Pi());
/*	
	*eta_emb=eta_rnd;
	*phi_emb=phi_rnd;
*/	

	if(jetType==0) //single-particle jet
	{	
		TLorentzVector v;
	   v.SetPtEtaPhiM(pT_emb, eta_rnd, phi_rnd, 0);
   	PseudoJet embeddedParticle=PseudoJet(v.Px(), v.Py(), v.Pz(), v.E());
	   embeddedParticle.set_user_index(99999);
		container.push_back(embeddedParticle);
	}
	else if(jetType==1) //pythia jet
	{
//cout<<"Memory used1: "<<memstat.Used()<<endl;
	
		bool charged=0; //we embed full jets since we need exact jet pT
  		TClonesArray *simarr = new TClonesArray("TParticle", 1000);
	
		//bool goodjet=0;
  		//while(!goodjet)
		//{
     	//simarr->Delete();
      Double_t phi_parton = phi_rnd;
      Double_t eta_parton = eta_rnd;
      Double_t theta = 2.0*TMath::ATan(TMath::Exp(-1*eta_parton));
      Double_t E = (1.03*pT_emb+0.5)/TMath::Sin(theta); //most of the jets (after constituent cuts) will have only 95% of the required pTemb => *1/0.95 = 1.05
      
		if(parton=="u")
   	   pythia->Py1ent(0, 2, E, theta, phi_parton); //u->jet
		else if (parton=="g")
	      pythia->Py1ent(0, 21, E, theta, phi_parton); //g->jet

      pythia->Pyexec();

      TLorentzVector partonlv(0, 0, 0, 0);
      partonlv.SetPtEtaPhiE(pT_emb, eta_parton, phi_parton, E);

      Int_t final = pythia->ImportParticles(simarr, "Final");
      Int_t nparticles = simarr->GetEntries();
     //cout<<"npart: "<<nparticles<<endl; 
      //TLorentzVector conelvPart; //particle level jet vector
  		TLorentzVector partlv; //particle level jet constituent vector

      Int_t goodparpart = 0;
		TLorentzVector foundlv(0., 0., 0., 0.);
      for(Int_t ipart = 0; ipart < nparticles; ipart++)
		{
			TParticle *particle = (TParticle*)simarr->At(ipart);

	  
	  		particle->Momentum(partlv);
	 		//cout<<ipart<<"| eta:"<<partlv.Eta()<<" phi: "<<partlv.Phi()<<endl; 
	  		partlv.SetPtEtaPhiM(partlv.Pt(), partlv.Eta(), partlv.Phi(), 0); //set M=0, probably not necessary

	  		Double_t pTpart = partlv.Pt();
    
 	  		if(pTpart < 0.2) continue; //undetectable soft particles
	 
	 		if(charged)
         {
          	Double_t charge = particle->GetPDG()->Charge();
          	if(!charge) continue;
         }

	  		Double_t eta = partlv.Eta();
		  	Double_t phi = partlv.Phi();
			Double_t M = 0;

	  		//conelvPart += partlv;
			PseudoJet embeddedParticle=PseudoJet(partlv.Px(), partlv.Py(), partlv.Pz(), partlv.E());
	   	embeddedParticle.set_user_index(99999);
			container.push_back(embeddedParticle);
			if(TMath::Abs(eta<1)) foundlv += partlv;

		}//particle loop [in particle level jet]
     	simarr->Delete();
		delete simarr;

		pT_emb=foundlv.Pt(); //update pT_emb value

	}//pythia jet
	return container;
}


//----------------------------------------------------------------------------- 
//Find embedded jet
//----------------------------------------------------------------------------- 
bool FindEmbeddedJet(vector<PseudoJet> constituents,float pT_emb)
{

		TLorentzVector foundlv(0., 0., 0., 0.);
    TLorentzVector constlv(0., 0., 0., 0.);
		int N=constituents.size();
    for(int iConst = 0; iConst < N; iConst++)
	   if(constituents[iConst].user_index() >= 90000) //embedded particles
      {
			float pT = constituents[iConst].perp();
			float eta = constituents[iConst].eta();
			float phi = constituents[iConst].phi();
			if(phi<0)phi=phi+2*TMath::Pi();
			float M = constituents[iConst].m();
			constlv.SetPtEtaPhiM(pT, eta, phi, M);
			foundlv += constlv;
		}

		if(foundlv.Pt() / pT_emb > 0.95) return true; //the jet is matched to the embedded one if it contains embedded particles which carry at least 90% energy of the embedded jet
		else return false;
}


//----------------------------------------------------------------------------- 
//Project track onto BEMC - similar to StEmcPosition::projTrack()
//----------------------------------------------------------------------------- 
Bool_t projectTrack(TVector3* atFinal, TVector3* momentumAtFinal, StPicoPhysicalHelix* helix, Double_t magField, Double_t radius = 229, Int_t option = 1)
{
    

     pair<double,double> pathLength = helix->pathLength(radius);

     Double_t s,s1,s2; 
     s=0;
     s1 = pathLength.first;
     s2 = pathLength.second;

     Bool_t goProj;
     goProj = kFALSE;
   
     if (finite(s1) == 0 && finite(s2) == 0) {cout << "track couldn't be projected" << endl; return kFALSE;} // Track couldn't be projected!
   
     if (option == 1)  // Selects positive path lenght to project track forwards along its helix relative to
                       // first point of track. The smaller solution is taken when both are positive
     {
       if (s1 >= 0 && s2 >= 0) {s = s1; goProj = kTRUE; }
       if (s1 >= 0 && s2 < 0) { s = s1; goProj = kTRUE; }
       if (s1 < 0 && s2 >= 0) { s = s2; goProj = kTRUE; }
     }
     
     if (option == -1) // Selects negative path lenght to project track backwards along its helix relative to
                       // first point of track. The smaller absolute solution is taken when both are negative 
     {
       if (s1 <= 0 && s2 <= 0) { s = s2; goProj = kTRUE; }
       if (s1 <= 0 && s2 > 0) { s = s1; goProj = kTRUE; }
       if (s1 > 0 && s2 <= 0) { s = s2; goProj = kTRUE; }
     }

     if (goProj) 
     {
       *atFinal = helix->at(s);
       *momentumAtFinal = helix->momentumAt(s, magField*tesla);
     }

 return goProj;
}

//----------------------------------------------------------------------------- 
//Find energy of neighboring towers - known eta, phi - obsolete
//----------------------------------------------------------------------------- 

Bool_t findNeighborTowers(StPicoDst* mPicoDst, StEmcGeom *mEmcGeom, double eta, double phi, int towerId, double *e1, double *e2/*, int *nextId1, int *nextId2*/)
{

	 //mEmcGeom = StEmcGeom::getEmcGeom("bemc");
   //int towerId = 0;
   int localTowerId = -1;
   //int localId1 = -1;
   //int localId2 = -1;
   double energy1 = 0, energy2 = 0;
   double energyTemp = 0;
   double dist1 = 1000, dist2 = 1000;
   double distTemp = 0;
   Float_t etaTemp = 0, phiTemp = 0;
 
   //if (mEmcGeom->getId(position.phi(), position.pseudoRapidity(), towerId) == 1) return kTRUE;
  
    for (int ieta = -1; ieta < 2; ++ieta) {
      for (int iphi = -1; iphi < 2; ++iphi) {
        localTowerId++;//loops from 0 to 8
				StEmcPosition *emcpos = new StEmcPosition();
        int nextTowerId = emcpos->getNextTowerId(towerId, ieta, iphi);
        if (nextTowerId < 1 || nextTowerId > 4800) continue;
		StPicoBTowHit *emcHit = mPicoDst->btowHit(nextTowerId-1);

        if (emcHit == 0) continue;
        if (emcHit->energy() < 0.2) continue; // don't include any noise tower


        if (ieta == 0 && iphi == 0) {
          mEmcGeom->getEta(nextTowerId, etaTemp);
          mEmcGeom->getPhi(nextTowerId, phiTemp);
          //ene[2] = emcHit->energy();
          //d[2] = position.pseudoRapidity() - etaTemp;
          //d[3] = position.phi() - phiTemp;
        }
        else {
          energyTemp = emcHit->energy();
          mEmcGeom->getEta(nextTowerId, etaTemp);
          mEmcGeom->getPhi(nextTowerId, phiTemp);
          distTemp = sqrt((etaTemp - eta) * (etaTemp - eta) + (phiTemp - phi) * (phiTemp - phi));

          // In case the new tower is closer to the matched tower
          // than the other closest one we swap them.
          // i.e. previously closest tower will become the second
          // closest tower
          if (distTemp < dist1) {
            dist2 = dist1;
            dist1 = distTemp;
            energy2 = energy1;
            energy1 = energyTemp;
            //localId2 = localId1;
            //localId1 = localTowerId;
          }
          else if (distTemp < dist2) {
            dist2 = distTemp;
            energy2 = energyTemp;
            //localId2 = localTowerId;
          }
        } //else
      } //for (int iphi = -1; iphi < 2; ++iphi)
    }
	*e1 = energy1;
	*e2 = energy2;

//	*nextId1 = localId1;
//	*nextId2 = localId2;

	return true;
}   

//----------------------------------------------------------------------------- 
//Find energy of neighboring towers - known position vector - obsolete
//----------------------------------------------------------------------------- 

Bool_t findNeighborTowers(StPicoDst* mPicoDst, StEmcGeom *mEmcGeom, TVector3 *position, int towerId, double *e1, double *e2/*, int *nextId1, int *nextId2*/)
{

	 //mEmcGeom = StEmcGeom::getEmcGeom("bemc");
   //int towerId = 0;
   int localTowerId = -1;
   //int localId1 = -1;
   //int localId2 = -1;
   double energy1 = 0, energy2 = 0;
   double energyTemp = 0;
   double dist1 = 1000, dist2 = 1000;
   double distTemp = 0;
   Float_t etaTemp = 0, phiTemp = 0;
 
   //if (mEmcGeom->getId(position.phi(), position.pseudoRapidity(), towerId) == 1) return kTRUE;
  
    for (int ieta = -1; ieta < 2; ++ieta) {
      for (int iphi = -1; iphi < 2; ++iphi) {
        localTowerId++;//loops from 0 to 8
				StEmcPosition *emcpos = new StEmcPosition();
        int nextTowerId = emcpos->getNextTowerId(towerId, ieta, iphi);
        if (nextTowerId < 1 || nextTowerId > 4800) continue;
		StPicoBTowHit *emcHit = mPicoDst->btowHit(nextTowerId-1);

        if (emcHit == 0) continue;
        if (emcHit->energy() < 0.2) continue; // don't include any noise tower


        if (ieta == 0 && iphi == 0) {
          mEmcGeom->getEta(nextTowerId, etaTemp);
          mEmcGeom->getPhi(nextTowerId, phiTemp);
          //ene[2] = emcHit->energy();
          //d[2] = position.pseudoRapidity() - etaTemp;
          //d[3] = position.phi() - phiTemp;
        }
        else {
          energyTemp = emcHit->energy();
          mEmcGeom->getEta(nextTowerId, etaTemp);
          mEmcGeom->getPhi(nextTowerId, phiTemp);
          distTemp = sqrt((etaTemp - position->PseudoRapidity()) * (etaTemp - position->PseudoRapidity()) + (phiTemp - position->Phi()) * (phiTemp - position->Phi()));

          // In case the new tower is closer to the matched tower
          // than the other closest one we swap them.
          // i.e. previously closest tower will become the second
          // closest tower
          if (distTemp < dist1) {
            dist2 = dist1;
            dist1 = distTemp;
            energy2 = energy1;
            energy1 = energyTemp;
            //localId2 = localId1;
            //localId1 = localTowerId;
          }
          else if (distTemp < dist2) {
            dist2 = distTemp;
            energy2 = energyTemp;
            //localId2 = localTowerId;
          }
        } //else
      } //for (int iphi = -1; iphi < 2; ++iphi)
    }
	*e1 = energy1;
	*e2 = energy2;

//	*nextId1 = localId1;
//	*nextId2 = localId2;

	return true;
}   



//----------------------------------------------------------------------------- 
//Find energy and position of a cluster surrounding track-matched tower 
//----------------------------------------------------------------------------- 

Bool_t findRealCluster(StPicoDst *mPicoDst, StEmcGeom *mEmcGeom, TVector3 *position, int towerId, double *ecluster, double *etacluster, double *phicluster, array<int, 9> *ids)
{
   int localTowerId = -1;
   double energy1 = 0, energy2 = 0;
   double energyTemp = 0;
   double dist1 = 1000, dist2 = 1000;
   double distTemp = 0;
   Float_t etaTemp = 0, phiTemp = 0;
		double cenergy = 0;
		array<double, 9> e;
		array<double, 9> phi;
		array<double, 9> eta;
		e.fill(0);
		phi.fill(0);
		eta.fill(0);

   //if (mEmcGeom->getId(position.phi(), position.pseudoRapidity(), towerId) == 1) return kTRUE;
  
    for (int ieta = -1; ieta < 2; ++ieta) {
      for (int iphi = -1; iphi < 2; ++iphi) {
        localTowerId++;//loops from 0 to 8
				StEmcPosition *emcpos = new StEmcPosition();
        int nextTowerId = emcpos->getNextTowerId(towerId, ieta, iphi);	
				bool bad = false;
				for (int i = 0; i<nBadTowers; i++) {if (BadTowerArr[i] == nextTowerId) bad = true;}
				if (bad) {
					//cout << "found bad tower " << nextTowerId <<" , throw cluster away" << endl; 
					for (int i=0; i<9; i++) ids->at(i) = -1; 
					return false;
					}
        if (nextTowerId < 1 || nextTowerId > 4800) continue;

				ids->at(localTowerId) = nextTowerId;
		StPicoBTowHit *emcHit = mPicoDst->btowHit(nextTowerId-1);

        if (emcHit == 0) continue;
        if (emcHit->energy() < 0.2) continue; // don't include any noise tower
			
				cenergy+=emcHit->energy();
				

 				if (ieta == 0 && iphi == 0) {
          mEmcGeom->getEta(nextTowerId, etaTemp);
          mEmcGeom->getPhi(nextTowerId, phiTemp);
          //ene[2] = emcHit->energy();
          //d[2] = position.pseudoRapidity() - etaTemp;
          //d[3] = position.phi() - phiTemp;
        }
        else {
          energyTemp = emcHit->energy();
          mEmcGeom->getEta(nextTowerId, etaTemp);
          mEmcGeom->getPhi(nextTowerId, phiTemp);
          distTemp = sqrt((etaTemp - position->PseudoRapidity()) * (etaTemp - position->PseudoRapidity()) + (phiTemp - position->Phi()) * (phiTemp - position->Phi()));

         // In case the new tower is closer to the matched tower
          // than the other closest one we swap them.
          // i.e. previously closest tower will become the second
          // closest tower
          if (distTemp < dist1) {
            dist2 = dist1;
            dist1 = distTemp;
            energy2 = energy1;
            energy1 = energyTemp;
            //localId2 = localId1;
            //localId1 = localTowerId;
          }
          else if (distTemp < dist2) {
            dist2 = distTemp;
            energy2 = energyTemp;
            //localId2 = localTowerId;
          }
        		} //else

			e[localTowerId] = emcHit->energy();
				//cout << towerId << " energy " << e[localTowerId] << " ieta " << ieta << " eta " << etaTemp << " iphi " << iphi << " phi " << phiTemp  << endl;
		//weighing phi and eta of individual towers by their energy
			phi[localTowerId] = phiTemp*e[localTowerId];
			eta[localTowerId] = etaTemp*e[localTowerId];
				//cout << towerId << " energy " << e[localTowerId] << " eta " << eta[localTowerId] << " phi " << phi[localTowerId] << endl;
      } //for (int iphi = -1; iphi < 2; ++iphi)
    }

		*ecluster = cenergy;
		if (cenergy > 0)
		{
		double sum = accumulate(phi.begin(), phi.end(), 0.0);

		*phicluster = sum/cenergy;
		//cout << accumulate(phi.begin(), phi.end(), sum) << " " << cenergy << endl;
		sum = accumulate(eta.begin(), eta.end(), 0.0);
		*etacluster = sum/cenergy;
		//cout << sum << " " << cenergy << endl;
		}

	e.fill(0);
	phi.fill(0);
	eta.fill(0);
	return true;
}   



//----------------------------------------------------------------------------- 
//find if tower in a cluster has been matched before
//----------------------------------------------------------------------------- 

bool doesExist(const vector< array<int,9> >&  v, int item, int *cluster){

    for (unsigned int row = 0; row < v.size(); row++) {
				//cout << row << endl;
       for (int col = 0; col < 9; col ++)
						{
						//cout << v[row][col] << " ";
							if (v[row][col] == item) 
							{//cout << "found tower: " << item << " in row: " << row << " column " << col << endl;
							//cout << v[row][col] << "";
							*cluster = row; 
            return true;
							}
						}
				//cout << endl;
    }

    return false;
}

//----------------------------------------------------------------------------- 
//find if matched tower is in a cluster
//----------------------------------------------------------------------------- 


bool doesExist(const vector< array<int,9> >&  v, int item){

    for (unsigned int row = 0; row < v.size(); row++) {
				//cout << row << endl;
       for (int col = 0; col < 9; col ++)
						{
						//cout << v[row][col] << " ";
							if (v[row][col] == item) 
					     return true;
					}
				//cout << endl;
    }

    return false;
}

//----------------------------------------------------------------------------- 
//find if tower exists in multiple clusters and save the cluster IDs (-1) into clusters vector
//----------------------------------------------------------------------------- 

bool doesExistInMultiple(const vector< vector<int> >&  Clusters, int item, vector<unsigned int> *clusters){
		bool found = false;
		
    for (unsigned int row = 0; row < Clusters.size(); row++) {
				//cout << row << endl;
       for (unsigned int col = 0; col < 9; col ++)
						{
							if (Clusters[row][col] == item) 
							{//cout << "found tower: " << item << " in cluster: " << row << " column " << col << endl;
							//cout << Clusters[row][col] << "";
							clusters->push_back(row); 
							found = true;
							}
						}
				//cout << endl;
    }
    return found;
}

//----------------------------------------------------------------------------- 
//Match track to BEMC cluster - maybe impossible from PicoDst - missing StEmcCollection
//----------------------------------------------------------------------------- 

/*StEmcCluster* StWiciED0PicoDstMaker::GetMatchingCluster(StPicoDst *mPicoDst, TVector3 *position, const StEmcGeom* geom, StEmcDetector* detector,
                                                        Double_t& selDist, Double_t& selDistPhi, Double_t& selDistZ, Int_t& nMatchCls) {
  //  Matches track and looks for the cluster that contanis the selected tower
  //  input:
  //   	position      - position of track after projection
  //    geom          - geometry of the detector
  //    detector      - detector object itself
  //  output:
  //    selDist       - distance between the track and the closest cluster on the detector surface
  //    selDistPhi    - distance in                         *etacluster = sum/energy;
                }
Phi plane between the track and the closest cluster on the detector surface
  //    selDistZ      - distance in Z direction between the track and the closest cluster on the detector surface
  //    nMatchCls     - number of clusters that have the tower which track matches to
  //  return:
  //    StEmcCluster* - Pointer to the closest StEmcCluster, but only if the projection went well and
  //                    the cluster with a tower that track matches to was found, otherwise NULL
  //StThreeVectorD        position, momentum;
  
  //projecting track to the surface of the given subdetector
  if(fEmcPosition->projTrack(&position, &momentum, track, fMagField, geom->Radius())==kFALSE) {
    cerr << "------> StWiciED0PicoDstMaker::GetMatchingCluster: Cannot project the track! ABORTING!!!" << endl;
    return NULL; //projection failed, aborting
  }
  TVector3 trkOnDet(position.x(), position.y(), position.z()); //conversion from STAR to ROOT format
  
  
  StEmcCluster* selCluster = NULL;
  TVector3 clusterPosition;
  Double_t     dist = 9999,    distPhi = 9999,    distZ = 9999;
  Double_t  minDist = 9999, minDistPhi = 9999, minDistZ = 9999;
  Int_t towMod = -1, towEta = -1, towSub = -1, towID = -1;
  //Int_t hitMod = -1, hitEta = -1, hitSub = -1;
  //Float_t phi, eta;
  
  Int_t assCl =  0;
  
  
  //finding a tower that track matches to (yup, tower, no matter what detector is in use here)
  if(geom->getBin(position->phi(), position->pseudoRapidity(), towMod, towEta, towSub) || towSub==-1) {
    //getBin returns 0 if eta is in the range of BEMC, 1 otherwise, an aditional test for sub is require though
    //cerr << "------> StWiciED0PicoDstMaker::GetMatchingCluster: Track out of BEMC! ABORTING!!!" << endl;
    return NULL; //when out of BEMC -> we don't want it!
  }

	StEmcCollection col* = new StEmcCollection();
  if(detector->detectorId()==kBarrelEmcTowerId) { //only when searching for tower cluster
    if(fBtowGeom->getId(towMod, towEta, towSub, towID)) { //getId returns 0 when everything's OK, 1 otherwise
      cerr << "------> StWiciED0PicoDstMaker::GetMatchingCluster: No SoftID! ABORTING!!!" << endl;
      return NULL; //when towID is broken --> Goodbye!
    }
    //checking if there is a tower with some signal in it
		StPicoBTowHit *emcHit = mPicoDst->btowHit(nextTowerId-1);
    if(fBtowHits[towID]==NULL) {
      cerr << "------> StWiciED0PicoDstMaker::GetMatchingCluster: Empty Tower! ABORTING!!!" << endl;
      return NULL; //if towID correct but no hit in tower --> Goodbye!
    }
  }
  
  //in every cluster...
  StSPtrVecEmcCluster&  clusters = detector->cluster()->clusters();
  for(StSPtrVecEmcClusterIterator cIt = clusters.begin(); cIt != clusters.end(); cIt++) {
    StPtrVecEmcRawHit& hits = (*cIt)->hit();
    //...checking hit by unless the match is right
    for(StPtrVecEmcRawHitIterator hIt = hits.begin(); hIt != hits.end(); hIt++) {
      //converting form detector m-e-s to BEMC tower m-e-s
      //if(detector->detectorId()!=kBarrelEmcTowerId) { //if det is not BEMC tower
      //  if(geom->getEta((*hIt)->module(), (Int_t)(*hIt)->eta(), eta) || geom->getPhi((*hIt)->module(), (Int_t)(*hIt)->sub(), phi))
      //    continue; //can't convert to eta, phi values (should not happen when the hit is correct)
      //  if(fBtowGeom->getBin(phi, eta, hitMod, hitEta, hitSub))
      //    continue; //can't convert to m-e-s of the BEMC tower (should not happen)
      //}
      //else { //for BEMC tower just asigning the values
      //  hitMod = (*hIt)->module();
      //  hitEta = (*hIt)->eta();
      //  hitSub = (*hIt)->sub();
      //}
      //if(towMod==hitMod && towEta==hitEta && towSub==hitSub) {
      //  assCl++;
        clusterPosition.SetPtEtaPhi(geom->Radius(), (*cIt)->eta(), (*cIt)->phi());
        GetDistance(geom->Radius(), position, clusterPosition, dist, distZ, distPhi);
        
        if(dist<minDist) { //checking if that's the cloasest cluster
          minDist = dist;
          minDistPhi = distPhi;
          minDistZ = distZ;
          selCluster = (*cIt);
          assCl = 1;
        }
        break;
      //} //if(hit match)
    } //for(hits)
  } //for(clusters)
  selDist    = minDist;
  selDistPhi = minDistPhi;
  selDistZ   = minDistZ;
  nMatchCls  = assCl;
  return selCluster;
} //end of GetMatchingCluster()
*/
//_______________________________________________________________________________________________
//_______________________________________________________________________________________________
//Find 2X2 cluster - obsolete
//_______________________________________________________________________________________________

int findClusterId(int towerId)
{
	int Rows = 1200;
	int Columns = 4;
	for(int row = 0; row < Rows; ++row)
	{
   for(int col = 0; col < Columns; ++col)
   {
      if(indexarr[row][col] == towerId) return row;
   }
	}
	return -1;
}

//_______________________________________________________________________________________________
//Find 2X2 cluster - obsolete
//_______________________________________________________________________________________________

Bool_t findCluster(StPicoDst *mPicoDst, StEmcGeom *mEmcGeom, int towerId, double *ecluster, double *etacluster, double *phicluster)
{
		double energy = 0;
		array<double, 4> e;
		array<double, 4> phi;
		array<double, 4> eta;
		e.fill(0);
		phi.fill(0);
		eta.fill(0);
		Float_t etaTemp = 0;
		Float_t phiTemp = 0;
		
		int cid = findClusterId(towerId);
		
		//look for the other towers in the cluster
		for (int i=0; i<4; i++)
		{
			int towid = indexarr[cid][i];
			bool bad = false;
			for (int i = 0; i < nBadTowers; i++) 
			{if (BadTowerArr[i] == towid) bad = true;}
			if (bad) {/*cout << "found bad tower " << towid <<" , throw cluster " << cid << "  away" << endl;*/ return false;}
			mEmcGeom->getEta(towid, etaTemp);
      mEmcGeom->getPhi(towid, phiTemp);


		//	if (towid == towerId) continue //no double counting of matched-tower energy
			StPicoBTowHit *emcHit = mPicoDst->btowHit(towid-1);
			if (emcHit == 0 || emcHit->energy() < 0.2) continue;
			
			e[i] = emcHit->energy();
			energy+=e[i];

		//weighing phi and eta of individual towers by their energy
			phi[i] = phiTemp*e[i];
			eta[i] = etaTemp*e[i];

		}

		*ecluster = energy;
		if (energy > 0)
		{
			double sum = accumulate(phi.begin(), phi.end(), 0.0);

			*phicluster = sum/energy;
			sum = accumulate(eta.begin(), eta.end(), 0.0);
			*etacluster = sum/energy;
		}

	e.fill(0);
	phi.fill(0);
	eta.fill(0);

		
		
	
	return true;
}
