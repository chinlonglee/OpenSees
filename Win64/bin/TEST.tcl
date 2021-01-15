wipe
model BasicBuilder -ndm 2 -ndf 2

node 1 0 0
node 2 1 0

# uniaxialMaterial Steel02 1 10 100000 .05 15 .925 .15
uniaxialMaterial Elastic 1 100

element truss 1 1 2 0.1 1

mass 2 10 10

fix 1 1 1
fix 2 0 1

timeSeries Path 1 -dt .02 -values {0.0 0.0 -0.1080 -0.1010 -0.0880 -0.0950 -0.1200 -0.1420 -0.1280 -0.1100 -0.0850 -0.0850 -0.1310 -0.1761 -0.1941 -0.1621 -0.1440 -0.1080 -0.0820 -0.0420 -0.0660 -0.1310 -0.1901 -0.1961 -0.0660 0.0300 0.1410 -0.0490 -0.1280 -0.1440 -0.2031 -0.2601 -0.3251 -0.3061 -0.1721 -0.1971 -0.1631 -0.1641 -0.0670 0.0250 0.1500 0.2361 0.2521 0.3361 0.4632 0.4922 0.4191 0.3591 0.2711 0.2351 0.3391 0.4121 0.5302 0.6392 0.7322 0.6522 0.5992 0.4001 0.4001 0.0630 -0.5152 -0.7873 -0.6032 -0.4842 -0.2501 -0.0590 0.1340 0.3081 0.4992 0.7102 0.9953 1.2204 1.5305 1.4505 1.1604 0.9353 0.8923 0.9263 0.8393 0.9013 0.9933 1.2104 0.3281 -1.4805 -2.0707 -1.9906 -2.0307 -1.8206 -1.7306 -1.7506 -1.7506 -1.8106 -1.6305 -1.3504 -1.0904 -0.7823 -0.4291 -0.0170 0.3601 0.7853 1.1604 1.6005 1.9606 2.4108 2.7309 3.0410 3.2010 3.4211 2.8209 2.3208 -1.2004 -2.3708 -1.6405 -1.8706 -1.1004 -0.7532 -0.1731 0.1130 0.5332 0.8953 1.1904 1.7606 0.5762 -2.6309 -1.5505 -1.7306 -1.0103 -0.5792 0.2371 -0.6702 -1.9806 -1.6405 -1.6906 -1.4805 -1.2304 -1.0003 -0.7512 -0.5232 -0.2711 -0.0440 0.1881 -0.0950 -0.4331 -0.8383 -0.9513 -0.7162 -0.5992 -0.3341 -0.1080 0.1851 0.4201 0.6732 -0.0970 -0.3721 -0.0400 0.0110 0.3441 0.5652 0.8833 1.1304 1.3604 0.2191 0.2411 0.6832 0.6892 1.3204 1.3504 2.0407 -0.9313 -1.3104 -0.6922 -0.5462 0.0720 0.6752 -1.0703 -1.4905 -1.0703 -1.1604 -0.7622 -0.5592 -0.2151 -0.1260 -0.6742 -0.3241 -0.3371 -0.1090 0.0170 0.2991 0.4882 0.6082 0.2221 -0.0320 -0.2451 0.0770 0.2111 0.5682 0.8263 1.2104 1.4805 1.7406 0.4211 0.0290 0.2591 0.2931 -0.0550 -0.1470 0.1430 0.2061 0.4992 0.6452 0.9573 1.1304 1.4505 1.6305 1.9506 1.8606 1.9806 1.7706 1.2504 -1.2104 -0.5422 -0.3841 -0.3111 -1.1204 -1.6605 -2.4608 -2.0307 -1.8406 -1.3204 -0.9603 -0.3251 0.1541 0.8163 1.3204 1.8206 -0.0580 -0.1691 0.2851 0.4471 0.9833 1.4205 1.8506 2.4608 1.6906 -1.3804 -0.9993 -1.0904 -0.9073 -0.4692 -1.2504 -2.1107 -1.6205 -1.6906 -1.3104 -1.1104 -0.7733 -0.5102 -0.5442 -1.2004 -1.2104 -1.1604 -1.1504 -0.7172 -0.5462 0.0640 -0.8043 -1.6305 -0.8593 -0.9613 -0.3961 -0.1470 0.3191 0.6482 0.8763 0.4722 0.1981 -0.0270 0.2921 0.4451 0.7853 1.0303 1.3504 1.6105 1.8606 1.2804 0.6402 0.2041 0.3141 0.3731 0.4962 0.2351 -0.0840 -0.1681 -0.1130 -0.2291 -0.2481 -0.1571 -0.0690 0.1470 0.3791 0.5792 0.2551 -0.0410 -0.4281 -0.1330 0.0950 0.2301 -0.1290 -0.0500 0.0800 0.2101 0.3801 0.5102 0.1571 -0.0320 -0.1110 0.0050 0.0760 0.0350 -0.0950 -0.0360 -0.0160 0.0380 0.0850 -0.0560 -0.3041 -0.4211 -0.2441 -0.2361 -0.1771 -0.1290 -0.0180 0.2031 -0.1080 -0.0910 -0.0340 -0.1060 -0.1110 -0.0990 -0.0020 0.0730 0.2351 0.3551 0.7052 0.7793 0.1841 -0.2631 -0.1240 -0.0420 0.1591 0.0480 -0.2191 -0.4672 -0.4281 -0.2161 -0.0430 0.1591 0.3201 0.4191 0.1230 -0.1601 -0.2041 -0.0820 -0.2061 -0.1370 -0.0550 0.0530 0.1340 0.2661 0.2321 0.0790 -0.0080 0.2001 0.4351 0.4922 0.1911 0.0920 -0.0220 -0.0210 0.0520 0.0930 0.2551 0.3681 0.5252 0.5412 0.4251 0.3981 0.5592 0.7562 0.3651 0.4111 0.0980 -0.2041 -0.2491 -0.4051 -0.4131 -0.4712 -0.4331 -0.4581 -0.0570 0.1781 -0.2081 -0.4922 -0.5302 -0.3621 -0.4051 -0.3081 -0.3161 -0.2651 -0.2651 -0.2691 -0.3451 -0.3091 -0.2171 -0.0780 0.0870 0.2811 0.3101 0.3581 0.3411 0.3581 0.2871 0.3051 0.1120 0.2141 0.1360 0.3841 -0.8613 -1.3504 -1.3404 -1.3504 -1.1904 -1.0403 -0.8293 -0.6512 -0.4441 -0.2581 -0.0600 -0.0910 -0.1821 -0.1470 0.0850 0.1631 0.0500 0.2641 0.5822 0.8673 1.2004 1.7006 1.1104 -1.1004 -0.3661 -0.4451 -0.2361 -0.9603 -0.6562 -0.5972 -0.6702 -0.5522 -0.0270 0.3781 1.0703 1.6705 0.9473 0.4081 0.6672 0.1320 -0.0950 -0.5202 -0.8273 -1.1504 -1.1504 -0.8033 -0.3691 0.0290 0.5452 1.1804 1.6105 -0.2701 0.0340 -0.0560 0.0200 0.1460 0.5372 0.7983 -0.2051 -0.5902 -0.1691 -0.1751 -0.0280 0.0740 0.3821 0.5672 0.7532 0.8013 0.5922 0.3041 0.0230 0.0640 -0.4061 -0.4511 -0.0790 0.1681 0.5672 0.0930 -0.0550 0.0440 -0.1230 -0.2821 -0.4371 -0.3521 -0.2551 -0.1110 0.2051 0.5192 0.8543 1.1404 0.7332 0.2371 -0.3681 -0.2711 -0.2171 -0.8733 -0.9733 -0.5892 -0.3361 0.0770 0.2591 0.5082 0.3611 0.0810 -0.0560 -0.2091 -0.3171 -0.2381 -0.3761 -0.5502 -0.7222 -0.8033 -0.5232 -0.3401 -0.0110 0.0650 -0.0370 -0.0050 -0.1681 -0.4101 -0.0800 0.0790 0.3741 0.6152 0.6652 0.2541 -0.0570 -0.4742 -0.3561 -0.2431 -0.0480 0.1260 0.3791 0.2411 -0.2271 -0.4281 -0.6792 -0.6612 -0.5902 -0.5132 -0.4081 -0.3091 -0.2661 -0.5412 -0.6282 -0.9083 -1.1104 -0.8813 -0.7703 -0.5822 -0.4732 -0.3331 -0.1991 0.0200 0.2111 0.4321 0.6132 0.7672 0.9333 1.0703 1.1304 1.1904 1.2504 1.3304 1.5905 1.8006 2.0407 1.2404 0.4421 -0.1400 -0.6662 -0.5552 -0.6932 -0.9843 -1.2504 -1.1804 -1.0503 -0.9203 -0.7432 -0.8093 -0.8503 -0.8603 -0.8633 -0.8733 -0.8683 -0.8853 -0.5372 0.0520 0.2151 0.2451 0.5802 0.3141 0.2361 0.4852 0.5892 0.5252 0.3551 0.1971 0.1991 0.4922 0.3431 0.2881 0.4321 0.2391 0.0880 0.0770 -0.1480 -0.0770 -0.0190 0.0750 0.0440 -0.1450 -0.3161 -0.2411 -0.0280 0.1821 0.4261 0.4391 0.5122 0.4662 0.4792 0.1931 0.2221 0.2741 0.3931 0.5042 0.5772 0.5882 0.8223 0.7973 0.9493 0.3451 0.0450 -0.1230 -0.3471 -0.4261 -0.4161 -0.2751 -0.2701 0.0740 0.4281 -0.2311 -0.3871 -0.0830 0.1390 0.4451 0.0270 -0.6972 -0.7963 -0.2511 -0.1350 0.0790 -0.1150 -0.2511 -0.3331 -0.2691 -0.3011 -0.2001 -0.0670 -0.0380 0.1050 0.2961 0.3441 0.9573 0.8983 0.1791 -0.3621 -0.9943 -0.8073 -0.7442 -0.5392 -0.3301 -0.1280 0.0310 0.1480 0.5082 -0.0220 -0.4892 -0.3581 -0.6912 -0.5162 -0.3711 0.0880 0.6322 0.8413 1.2804 1.3905 1.1904 0.7512 0.2251 -0.0880 -0.2271 0.0740 0.1811 0.5442 0.3991 0.0450 -0.0820 -0.1851 -0.0200 0.0060 -0.1170 -0.2101 -0.3031 -0.5122 -0.7272 -0.5792 -0.2661 -0.1781 0.0400 0.0980 0.1370 0.2211 0.4371 0.0910 -0.5482 -0.5552 -0.2431 -0.0810 0.2501 0.4101 0.1821 -0.0270 -0.2431 -0.0150 0.2471 0.4822 0.7833 0.6222 0.3311 -0.0140 -0.1951 -0.2471 -0.2121 -0.1100 0.0500 0.2411 -0.0340 -0.2161 -0.4712 -0.3631 -0.1951 -0.0180 0.1701 -0.0800 0.0050 0.2301 0.3741 0.6012 0.5162 0.4321 0.3441 0.5052 0.6532 0.6832 0.1721 -0.1701 -0.5272 -0.6642 -0.3871 -0.2221 -0.0330 0.1190 -0.1280 -0.3511 -0.5142 -0.3351 -0.2181 -0.0120 0.1420 0.0700 -0.0630 -0.1200 -0.3221 -0.3461 -0.0910 0.0730 0.3091 0.4722 0.6032 0.5762 0.3301 -0.0730 -0.7773 -0.6082 -0.4381 -0.2091 0.0310 0.3501 0.2931 0.1210 0.3381 0.3171 0.2541 0.2061 0.1981 0.1741 0.0210 -0.1440 -0.3431 -0.3391 -0.1450 -0.0280 0.1701 -0.0960 -0.2551 -0.2791 -0.3881 -0.2421 -0.2151 -0.1821 -0.1741 -0.0380 -0.0270 -0.1851 -0.1230 0.0870 0.3431 0.6952 0.9103 0.8533 0.7602 0.5132 0.1861 0.0150 -0.1901 -0.1510 -0.0730 0.0210 0.1290 0.2151 0.0240 -0.1240 -0.3291 -0.5192 -0.7082 -0.5792 -0.4622 -0.3071 -0.1450 -0.0090 -0.1801 -0.3181 -0.4652 -0.3911 -0.3451 -0.3161 -0.4351 -0.4912 -0.4752 -0.4201 -0.3611 -0.2771 -0.2581 -0.1390 -0.0680 0.5072 0.7222 0.8783 0.7823 0.7652 0.4391 0.0800 0.0130 -0.1260 -0.0150 0.0300 0.1040 0.1040 0.1931 0.2051 0.0740 -0.0560 -0.0720 0.0700 0.1060 0.1470 -0.0090 -0.1591 -0.1871 -0.0070 0.1551 0.1050 -0.1150 -0.3021 -0.3091 -0.0950 -0.0580 0.0040 0.0200 0.0500 0.0570 0.0970 0.1340 0.1771 0.2181 0.2611 0.3021 0.3461 0.3861 0.4742 0.3931 0.2381 0.1150 -0.0790 -0.1240 0.0540 0.0270 -0.2501 -0.5662 -0.6302 -0.5912 -0.4131 -0.0680 0.2721 0.2771 -0.0210 -0.0600 -0.1100 -0.2211 -0.4161 -0.5192 -0.2221 0.0300 0.0790 0.1390 0.1711 0.2531 0.3231 0.3911 0.1641 -0.1360 -0.3231 -0.2911 -0.2871 -0.3041 -0.3391 -0.2451 -0.0760 0.1250 0.3761 0.4021 0.2451 0.1561 -0.0400 -0.1530 -0.2891 -0.3161 -0.1110 0.0940 0.3351 0.5762 0.4241 0.1430 -0.0070 -0.1350 -0.2701 -0.3411 -0.3571 -0.3961 -0.4021 -0.4882 -0.4802 -0.4061 -0.4071 -0.3511 -0.1871 -0.0570 0.0440 -0.0190 -0.0720 -0.1691 -0.1150 0.1260 0.3581 0.6542 0.7162 0.7622 0.7392 0.6282 0.4842 0.2641 -0.0440 -0.2881 -0.3841 -0.4922 -0.4281 -0.4161 -0.2761 -0.0520 0.2371 0.4261 0.6042 0.4521 0.2841 0.1260 -0.0540 -0.2751 -0.4231 -0.1751 0.0010 0.2541 0.4812 0.6412 0.5542 0.4131 0.1841 -0.0480 -0.3031 -0.5312 -0.7082 -0.9283 -0.8633 -0.6312 -0.3761 0.0870 0.3091 0.5892 0.6142 0.3851 0.3511 0.3111 0.2461 0.0190 -0.1981 -0.1581 -0.0150 0.0990 0.2861 0.4081 0.5632 0.5312 0.3141 0.1651 -0.0240 -0.1891 -0.2761 -0.3711 -0.4501 -0.5342 -0.4832 -0.3791 -0.2961 -0.1961 -0.1841 -0.1591 -0.1190 -0.0530 0.0020 0.0590 -0.0230 -0.1120 -0.2051 -0.3221 -0.3881 -0.3411 -0.2871 -0.3281 -0.4071 -0.4872 -0.5632 -0.6442 -0.5552 -0.4461 -0.0080 0.2531 0.4111 0.6442 0.5792 0.4742 0.3841 0.3851 0.3401 0.3571 0.0080 -0.2541 -0.4601 -0.4712 -0.2221 -0.0650 0.1641 0.3551 0.5042 0.3711 0.2801 0.1581 0.0410 -0.0270 0.0210 0.0440 0.0990 0.1350 0.0940 0.0580 0.0220 -0.0310 -0.0310 0.0040 -0.0250 -0.1240 -0.2351 -0.4061 -0.5302 -0.7012 -0.3231 -0.0520 0.0940 0.3281 0.4782 0.5092 0.3581 0.3421 0.3061 0.2851 0.2631 0.1681 -0.0050 -0.2121 -0.3651 -0.3101 -0.2971 -0.2801 -0.2371 -0.2661 -0.3081 -0.3661 -0.3571 -0.3081 -0.1931 0.0190 0.1961 0.1601 0.1290 0.1400 0.1100 0.1080 0.0920 0.0890 0.0190 -0.1310 -0.2471 -0.4361 -0.4321 -0.3001 -0.1921 -0.0480 0.0970 0.1681 0.1480 0.1731 0.0780 -0.0580 -0.2151 -0.2341 0.0600 0.2621 0.2691 0.0840 -0.0410 -0.2271 -0.0760 0.0310 0.1821 0.3101 0.4792 0.4591 0.1661 -0.0580 -0.3961 -0.4441 -0.2411 -0.1010 0.2171 0.2611 0.1210 0.0080 -0.1681 -0.3821 -0.5662 -0.7813 -0.6192 -0.2561 0.0440 0.4561 0.7833 1.1004 0.9533 0.4892 0.1220 -0.3871 -0.8463 -1.2304 -0.8643 -0.6162 -0.2281 0.1500 0.4772 0.3231 0.2481 0.0950 0.2271 0.3381 0.5052 0.5942 0.5192 0.5522 0.5952 0.6172 0.5932 0.4662 0.1581 -0.0770 -0.2301 -0.3031 -0.3291 -0.3641 -0.4902 -0.6602 -0.7232 -0.7843 -0.8463 -0.5662 -0.2091 0.1280 0.3811 0.4882 0.3631 0.2501 0.1581 0.3771 0.5942 0.6302 0.3701 0.2171 0.0400 0.0230 -0.0620 -0.2861 -0.3801 -0.4571 -0.5052 -0.7773 -0.6282 -0.1891 0.1891 0.5092 0.7442 0.7492 0.5622 0.4551 0.5252 0.6352 0.7312 0.8513 0.9273 0.9753 0.9293 0.8143 0.2641 -0.2891 -0.9013 -1.3604 -0.8593 -0.4682 0.0810 0.5942 0.8233 0.5562 0.3011 -0.0110 -0.3801 -0.3471 -0.0200 0.1801 0.6332 1.0103 1.1204 1.1104 0.9403 0.8053 0.6082 0.5082 0.2061 -0.2011 -0.5722 -1.0203 -1.3304 -1.2304 -1.2104 -1.0103 -0.7622 -0.5482 -0.4892 -0.4051 -0.2911 -0.1761 -0.1410 -0.0980 -0.0500 0.0190 0.0800 0.0350 -0.0290 -0.0610 -0.0340 -0.0130 0.0760 0.2051 0.3271 0.4501 0.5772 0.5622 0.4642 0.3001 0.1070 -0.0040 -0.0330 -0.0690 -0.1010 -0.2011 -0.1881 0.0250 0.1340 0.2431 0.2661 0.2651 0.1551 -0.0190 -0.0950 -0.2161 -0.1671 -0.1200 -0.0640 -0.1290 -0.1631 -0.1931 -0.2421 -0.2361 -0.1751 -0.1240 -0.1851 -0.2651 -0.3231 -0.3361 -0.4541 -0.4301 -0.3341 -0.2131 -0.0690 0.0300 0.0030 -0.0930 -0.0890 -0.1500 -0.1641 -0.2381 -0.3231 -0.4211 -0.4571 -0.3971 -0.3491 -0.2581 -0.1721 -0.0200 0.1561 0.2841 0.3621 0.3541 0.2691 0.1010 -0.0450 -0.1250 -0.2451 -0.2291 -0.1260 -0.0680 0.0180 0.0930 0.2001 0.2861 0.3651 0.3111 0.1811 0.0240 -0.1561 -0.3191 -0.2191 -0.1180 0.0150 0.1530 0.2981 0.2431 0.1360 0.1020 0.0240 -0.0100 -0.0230 -0.0360 -0.0480 -0.0890 -0.1551 -0.1480 -0.0760 0.0000 -0.0200 -0.1480 -0.2251 -0.3741 -0.3651 -0.2511 -0.1641 -0.0140 0.1530 0.3021 0.3931 0.4061 0.3851 0.3291 0.1831 0.1250 0.0700 -0.0010 -0.0680 -0.1300 -0.1280 -0.1150 -0.1020 -0.0840 -0.1410 -0.2051 -0.2691 -0.3501 -0.3591 -0.2971 -0.2271 0.0180 0.0990 0.2031 0.2561 0.1991 0.1601 0.1070 0.1561 0.1981 0.2421 0.1410 0.0550 -0.0660 -0.0970 -0.0260 0.0160 0.1000 0.0940 0.0390 0.0000 -0.0660 -0.0990 -0.1090 -0.1150 -0.1310 -0.1681 -0.1991 -0.1040 -0.0100 0.0550 0.1290 0.2001 0.1991 0.1691 0.1430 0.1090 0.1280 0.1520 0.1370 0.1170 0.0970 0.0600 -0.0280 -0.0500 -0.0840 -0.1781 -0.3211 -0.4201 -0.5182 -0.4722 -0.3941 -0.2901 -0.1100 0.0500 0.1320 0.1571 0.1771 0.1951 0.2551 0.3311 0.3281 0.2521 0.1731 0.0450 -0.0850 -0.1901 -0.2291 -0.3041 -0.2771 -0.2271 -0.1741 -0.1210 -0.1260 -0.1290 -0.0800 -0.0260 0.0390 0.0650 0.0550 0.0500 0.0640 0.1260 0.1791 0.2431 0.3071 0.2981 0.2511 0.2161 0.1631 0.1921 0.2341 0.2821 0.3011 0.1961 0.1060 -0.0390 -0.1631 -0.3221 -0.3351 -0.2191 -0.1480 -0.0100 0.0110 -0.0530 -0.1010 -0.1260 -0.1430 -0.1290 -0.1040 -0.0710 -0.0180 0.0330 0.0850 0.1581 0.2391 0.3191 0.3411 0.3181 0.2131 0.0730 0.0090 -0.0150 -0.0470 -0.0750 -0.1210 -0.1561 -0.1160 -0.0560 -0.0060 -0.0030 -0.0030 0.0140 0.0460 0.0720 0.0860 0.0980 0.1240 0.1470 0.1721 0.2001 0.2561 0.3171 0.2871 0.2311 0.1050 -0.0010 -0.0110 -0.0360 -0.0530 -0.0830 -0.0520 -0.0070 0.0370 0.0960 0.1551 0.2051 0.1430 0.0730 0.0210 -0.0220 -0.0700 -0.1000 -0.0750 -0.0540 -0.0290 -0.0150 0.0000 -0.0010 -0.0060 -0.0110 -0.0100 -0.0030 0.0010 0.0160 0.0530 0.0860 0.1260 0.1541 0.1310 0.1020 0.0560 0.0060 -0.0400 -0.0980 -0.0960 -0.0460 -0.0070 0.0310 0.0450 0.0680 0.0520 0.0450 0.0180 -0.0020 -0.0290 -0.0280 -0.0170 -0.0060 -0.0050 -0.0170 -0.0220 -0.0360 -0.0320 -0.0070 0.0140 0.0410 0.0650 0.0630 0.0520 0.0060 -0.0450 -0.0760 -0.0640 -0.0650 -0.1070 -0.1611 -0.1741 -0.1170 -0.0700 -0.0360 -0.0250 -0.0020 0.0040 0.0360 0.0820 0.1150 0.0700 0.0380 -0.0100 -0.0250 -0.0340 -0.0410 -0.0450 -0.0650 -0.1200 -0.1170 -0.1080 -0.0980 -0.0810 -0.0760 -0.1100 -0.1380 -0.1661 -0.1440 -0.1300 -0.1110 -0.1080 -0.0980 -0.0980 -0.0630 -0.0020 0.0540 0.1100 0.1520 0.1971 0.1701 0.1400 0.1030 0.0640 0.0100 -0.0670 -0.1050 -0.1220 -0.1380 -0.1661 -0.2151 -0.2691 -0.2741 -0.2171 -0.1781 -0.1160 -0.0610 0.0000 0.0540 0.0990 0.1440 0.1871 0.2321 0.2621 0.2551 0.2271 0.1711 0.1270 0.0540 -0.0470 -0.1490 -0.2101 -0.2541 -0.2521 -0.2451 -0.2261 -0.1791 -0.1360 -0.0820 -0.0140 0.0440 0.0890 0.0970 0.1000 0.1000 0.1040 0.1360 0.1210 0.1000 0.0720 0.0190 -0.0280 -0.0830 -0.1310 -0.1901 -0.1941 -0.1390 -0.1000 -0.0440 -0.0440 -0.0400 -0.0450 -0.0490 -0.0390 -0.0210 -0.0220 -0.0160 -0.0160 -0.0120 -0.0120 0.0280 0.0680 0.1150 0.0960 0.0560 0.0170 -0.0310 -0.0780 -0.0920 -0.0920 -0.0950 -0.0920 -0.0860 -0.0520 -0.0220 0.0140 0.0380 0.0340 0.0340 0.0290 0.0400 0.0540 0.0670 0.0850 0.0870 0.0770 0.0700 0.0610 0.0750 0.0870 0.1020 0.1110 0.0900 0.0750 0.0530 0.0620 0.0790 0.0970 0.1180 0.1350 0.1050 0.0770 0.0440 0.0120 -0.0210 -0.0420 -0.0680 -0.0870 -0.1200 -0.1591 -0.1981 -0.2331 -0.2581 -0.2701 -0.2781 -0.2871 -0.2781 -0.2701 -0.2581 -0.2511 -0.2191 -0.1731 -0.1290 -0.1090 -0.1050 -0.0950 -0.0910 -0.0720 -0.0590 -0.0400 -0.0230 0.0180 0.0350 0.0420 0.0480 0.0430 0.0430 0.0380 0.0360 0.0320 0.0470 0.0630 0.0850 0.1100 0.1340 0.1601 0.1811 0.1731 0.1681 0.1591 0.1470 0.1080 0.0680 0.0320 0.0090 0.0190 0.0280 0.0400 0.0540 0.0690 0.0950 0.1170 0.1430 0.1430 0.1250 0.1220 0.1240 0.1260 0.1300 0.1330 0.1060 0.0790 0.0460 0.0200 0.0050 -0.0110 -0.0250 -0.0110 0.0180 0.0430 0.0760 0.1030 0.1330 0.1010 0.0620 0.0180 -0.0270 -0.0800 -0.0890 -0.0540 -0.0300 0.0060 0.0120 0.0220 0.0270 0.0320 0.0480 0.0690 0.0890 0.1120 0.1330 0.1400 0.1370 0.1400 0.1350 0.1390 0.1020 0.0300 -0.0320 -0.0630 -0.0850 -0.1120 -0.1571 -0.1971 -0.1881 -0.1821 -0.1701 -0.1561 -0.1420 -0.1290 -0.1250 -0.1200 -0.1160 -0.1140 -0.1060 -0.0640 -0.0210 0.0020 -0.0040 -0.0040 -0.0140 -0.0190 -0.0080 0.0110 0.0270 0.0500 0.0530 0.0420 0.0300 0.0110 -0.0050 -0.0240 -0.0200 -0.0070 0.0040 0.0190 0.0330 0.0480 0.0570 0.0650 0.0740 0.0810 0.0910 0.1270 0.1691 0.1440 0.1050 0.0700 0.0350 0.0020 -0.0320 -0.0640 -0.0950 -0.0910 -0.0830 -0.0730 -0.0510 -0.0300 -0.0040 0.0330 0.0520 0.0440 0.0400 0.0300 0.0200 0.0100 0.0010 -0.0100 -0.0140 0.0020 0.0160 0.0300 0.0230 0.0160 0.0160 0.0260 0.0320 0.0450 0.0460 0.0210 -0.0030 -0.0300 -0.0600 -0.0870 -0.1160 -0.1070 -0.0940 -0.0780 -0.0680 -0.0630 -0.0570 -0.0500 -0.0410 -0.0330 -0.0240 -0.0140 -0.0050 0.0040 0.0130 0.0210 0.0010 -0.0180 -0.0400 -0.0640 -0.0870 -0.1100 -0.1290 -0.1510 -0.1651 -0.1741 -0.1841 -0.1941 -0.2021 -0.2131 -0.1981 -0.1731 -0.1350 -0.0650 -0.0140 0.0200 0.0550 0.0820 0.0890 0.1020 0.0930 0.0800 0.0610 0.0330 0.0030 -0.0240 -0.0570 -0.0450 -0.0130 0.0130 0.0330 0.0480 0.0540 0.0580 0.0620 0.0650 0.0690 0.0720 0.0750 0.0770 0.0730 0.0690 0.0640 0.0580 0.0540 0.0480 0.0610 0.0790 0.0960 0.1140 0.1260 0.1390 0.1310 0.1230 0.1140 0.1030 0.0870 0.0720 0.0560 0.0400 0.0240 0.0070 0.0040 0.0080 0.0090 0.0230 0.0420 0.0520 0.0610 0.0720 0.0790 0.0890 0.0950 0.0950 0.0930 0.0890 0.0250 0.0000 -0.0050 -0.0200 -0.0240 -0.0330 -0.0480 -0.0640 -0.0800 -0.0970 -0.1140 -0.1280 -0.1060 -0.0890 -0.0660 -0.0430 -0.0190 -0.0110 -0.0090 -0.0050 -0.0040 0.0050 0.0140 0.0240 0.0350 0.0430 0.0440 0.0390 0.0130 -0.0080 -0.0360 -0.0470 -0.0550 -0.0620 -0.0680 -0.0730 -0.0800 -0.0750 -0.0640 -0.0280 0.0050 0.0460 0.0650 0.0600 0.0600 0.0540 0.0490 0.0440 0.0390 0.0290 0.0000 -0.0260 -0.0520 -0.0420 -0.0370 -0.0250 -0.0150 -0.0020 0.0020 0.0000 -0.0010 -0.0040 -0.0070 -0.0060 -0.0010 0.0030 0.0070 -0.0010 -0.0070 -0.0150 -0.0240 -0.0290 -0.0280 -0.0280 -0.0250 -0.0240 -0.0210 -0.0270 -0.0440 -0.0580 -0.0760 -0.0840 -0.0820 -0.0830 -0.0750 -0.0500 -0.0280 -0.0020 0.0180 0.0180 0.0210 0.0210 0.0210 0.0310 0.0430 0.0490 0.0490 0.0490 0.0480 0.0470 0.0520 0.0670 0.0790 0.0950 0.0970 0.0960 0.0960 0.0930 0.0910 0.0890 0.0870 0.0840 0.0820 0.0800 0.0770 0.0620 0.0470 0.0350 0.0360 0.0340 0.0350 0.0360 0.0320 0.0190 0.0150 0.0140 0.0120 0.0120 0.0110 0.0000 -0.0100 -0.0200 -0.0250 -0.0290 -0.0330 -0.0370 -0.0480 -0.0630 -0.0770 -0.0920 -0.1070 -0.1220 -0.1380 -0.1430 -0.1200 -0.1030 -0.0780 -0.0670 -0.0630 -0.0580 -0.0560 -0.0530 -0.0500 -0.0480 -0.0450 -0.0420 -0.0390 -0.0310 -0.0260 -0.0110 0.0100 0.0310 0.0460 0.0480 0.0540 0.0550 0.0630 0.0840 0.1010 0.0840 0.0710 0.0530 0.0330 0.0280 0.0290 0.0280 0.0270 0.0120 -0.0010 -0.0100 -0.0110 -0.0160 -0.0150 -0.0310 -0.0680 -0.0990 -0.1020 -0.1050 -0.1060 -0.0990 -0.0890 -0.0790 -0.0740 -0.0720 -0.0680 -0.0670 -0.0630 -0.0400 -0.0130 0.0130 0.0430 0.0710 0.0990 0.0860 0.0630 0.0410 0.0130 -0.0140 -0.0350 -0.0370 -0.0440 -0.0460 -0.0490 -0.0510 -0.0540 -0.0560 -0.0580 -0.0460 -0.0290 -0.0130 0.0070 0.0170 0.0190 0.0250 0.0250 0.0350 0.0500 0.0640 0.0800 0.0940 0.0970 0.1020 0.1050 0.1080 0.1100 0.1040 0.0990 0.0930 0.0860 0.0790 0.0720 0.0630 0.0550 0.0460 0.0380 0.0310 0.0270 0.0210 0.0170 0.0120 0.0070 0.0030 0.0080 0.0130 0.0190 0.0220 0.0090 -0.0030 -0.0170 -0.0310 -0.0470 -0.0530 -0.0480 -0.0430 -0.0380 -0.0250 -0.0090 0.0070 0.0170 0.0220 0.0350 0.0610 0.0830 0.1120 0.1070 0.0870 0.0710 0.0470 0.0270 0.0040 -0.0160 -0.0240 -0.0330 -0.0460 -0.0690 -0.0840 -0.0880 -0.0930 -0.0950 -0.0970 -0.0990 -0.1010 -0.1020 -0.0980 -0.0950 -0.0900 -0.0870 -0.0890 -0.0910 -0.0930 -0.0960 -0.0970 -0.0880 -0.0810 -0.0710 -0.0680 -0.0680 -0.0670 -0.0670 -0.0660 -0.0660 -0.0640 -0.0560 -0.0500 -0.0430 -0.0340 -0.0230 -0.0130 -0.0020 0.0100 0.0170 0.0090 0.0030 -0.0050 -0.0140 -0.0070 0.0060 0.0160 0.0130 0.0120 0.0080 0.0050 0.0000 0.0020 0.0090 0.0140 0.0210 0.0220 0.0100 0.0000 -0.0120 -0.0240 -0.0240 -0.0190 -0.0160 -0.0180 -0.0220 -0.0270 -0.0300 -0.0260 -0.0220 -0.0220 -0.0320 -0.0410 -0.0510 -0.0460 -0.0430 -0.0380 -0.0320 -0.0260 -0.0200 -0.0190 -0.0180 -0.0150 0.0130 0.0410 0.0670 0.0700 0.0760 0.0780 0.0800 0.0810 0.0840 0.0780 0.0610 0.0480 0.0310 0.0290 0.0300 0.0320 0.0350 0.0370 0.0240 0.0070 -0.0080 -0.0280 -0.0350 -0.0220 -0.0130 0.0010 0.0030 0.0010 -0.0010 -0.0050 -0.0080 -0.0110 -0.0150 -0.0120 -0.0080 -0.0040 0.0010 0.0050 0.0080 0.0080 0.0090 0.0100 0.0220 0.0330 0.0460 0.0590 0.0720 0.0790 0.0620 0.0500 0.0320 0.0160 -0.0040 0.0000 0.0300 0.0520 0.0870 0.0820 0.0310 -0.0020 -0.0620 -0.0280 0.0580 0.1300 0.2211 0.1851 0.0940 0.0910 0.0770 0.0900 0.0160 -0.0840 -0.2001 -0.1480 -0.0070 0.1070 0.0990 0.0480 0.0030 -0.0690 -0.1070 -0.0150 0.0690 0.1380 0.0670 0.0200 0.0350 0.0900 0.0860 0.0520 0.0240 0.0160 0.0130 -0.0110 -0.0340 -0.0670 -0.0540 -0.0200 0.0140 0.0480 0.0970 0.0580 -0.0480 -0.1190 -0.0950 -0.0600 -0.0240 0.0320 0.0380 0.0140 -0.0010 -0.0070 0.0060 0.0100 0.0240 0.0130 0.0040 -0.0080 -0.0270 -0.0550 -0.0840 -0.0930 -0.0610 -0.0360 -0.0020 0.0330 0.0370 -0.0140 -0.0450 -0.1030 -0.0770 -0.0230 0.0310 0.0410 0.0260 0.0130 0.0150 0.0610 0.0730 0.0740 0.0760 0.0730 0.0710 0.0570 0.0380 0.0210 -0.0020 -0.0330 -0.0640 -0.0680 -0.0510 -0.0540 -0.0540 -0.0570 -0.0330 0.0000 0.0330 0.0700 0.0960 0.0850 0.0800 0.0680 0.0750 0.0770 0.0850 0.0610 0.0370 0.0070 -0.0010 0.0040 0.0050 -0.0100 -0.0290 -0.0460 -0.0690 -0.0650 -0.0440 -0.0280 -0.0070 -0.0230 -0.0330 -0.0510 -0.0600 -0.0560 -0.0540 -0.0500 -0.0660 -0.0910 -0.1140 -0.1400 -0.1380 -0.1400 -0.1210 -0.0960 -0.0740 -0.0460 -0.0340 -0.0420 -0.0450 -0.0550 -0.0630 -0.0720 -0.0800 -0.0730 -0.0520 -0.0370 -0.0260 -0.0140 0.0000 0.0000}

pattern UniformExcitation 1 1 -accel 1  -fact 10

constraints Plain
numberer Plain
algorithm Newton
system LeeSparse
integrator LeeNewmarkFullKC 0.5 0.25 -type0 1 .05 -type1 10 .03 3
analysis Transient

#analysis Transient
test NormDispIncr 1E-12 20 1

analyze 100 0.02