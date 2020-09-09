function [model, data] = GenerateModel_JeCh(data)
% version 1

% File extracted from JeCh_200.s2mMod using Spacial_v2_Modelizer.py

% Informations générales
% root_actuated	0
% external_forces	0

import casadi.*

model.name = 'JeCh_200';

model.NB = 37;

model.bodyN_name = {'Tx_pelvis','Ty_pelvis','Tz_pelvis','Rx_pelvis','Ry_pelvis','Rz_pelvis','Rx_thorax','Ry_thorax','Rz_thorax','Rx_tete','Ry_tete','Rz_tete','Ry_epauled','Rz_epauled','Rx_brasd','Ry_brasd','Rz_brasd','Rx_abrasd','Rz_abrasd','Rx_maind','Ry_maind','Ry_epauleg','Rz_epauleg','Rx_brasg','Ry_brasg','Rz_brasg','Rx_abrasg','Rz_abrasg','Rx_maing','Ry_maing','Rx_cuissed','Ry_cuissed','Rz_cuissed','Rx_jambed','Rx_piedd','Rz_piedd','Rx_cuisseg','Ry_cuisseg','Rz_cuisseg','Rx_jambeg','Rx_piedg','Rz_piedg'};
model.jtype = {'R','Rx','Ry','Rz','Rx','Ry','Rz','Ry','Rz','Rx','Ry','Rz','Rx','Rz','Rx','Ry','Ry','Rz','Rx','Ry','Rz','Rx','Rz','Rx','Ry','Rx','Ry','Rz','Rx','Rx','Rz','Rx','Ry','Rz','Rx','Rx','Rz'};
model.parent = [0,1,2,3,4,5,6,4,8,9,10,11,12,13,14,15,4,17,18,19,20,21,22,23,24,1,26,27,28,29,30,1,32,33,34,35,36];

model.Xtree = {eye(6),inv(pluho([0.9998112728 -0.0152768896 0.0120014726 -0.0105445918;0.0159703143 0.9980717122 -0.0599817167 -0.0594876700;-0.0110619962 0.0601620638 0.9981273257 0.2160444777;0.0000000000 0.0000000000 0.0000000000 1.0000000000])),eye(6),eye(6),inv(pluho([0.9946283617 -0.0053860420 0.1033702688 -0.0014793674;0.0915078209 0.5125237080 -0.8537832087 -0.0066127196;-0.0483812012 0.8586561822 0.5102634810 0.2994318521;0.0000000000 0.0000000000 0.0000000000 1.0000000000])),eye(6),eye(6),inv(pluho([0.9926071942 0.0960671994 0.0741758111 0.0323399606;-0.0905310952 0.9930877165 -0.0747054765 0.0167688444;-0.0808398328 0.0674379760 0.9944430807 0.1699134093;0.0000000000 0.0000000000 0.0000000000 1.0000000000])),eye(6),inv(pluho([0.7407461168 -0.6039532562 -0.2941694317 0.1612546624;0.6178599746 0.7843909783 -0.0545879575 -0.0000019510;0.2637124230 -0.1413197002 0.9541931147 0.0000024963;0.0000000000 0.0000000000 0.0000000000 1.0000000000])),eye(6),eye(6),inv(pluho([0.1785044571 -0.9826780529 -0.0497996259 0.0000006432;0.8537352156 0.1798460622 -0.4886630493 0.0000007417;0.4891547203 0.0447128379 0.8710501832 -0.2929873739;0.0000000000 0.0000000000 0.0000000000 1.0000000000])),eye(6),inv(pluho([0.9779455590 0.0350420881 0.2058993337 0.0000003366;0.0484240070 0.9209215480 -0.3867280926 0.0000005075;-0.2031688933 0.3881694915 0.8989142601 -0.2540782785;0.0000000000 0.0000000000 0.0000000000 1.0000000000])),eye(6),inv(pluho([0.9840921936 -0.1440107224 -0.1040358888 -0.0284824291;0.1366368693 0.9877881205 -0.0748665146 0.0276167488;0.1135469960 0.0594604144 0.9917517526 0.1667404889;0.0000000000 0.0000000000 0.0000000000 1.0000000000])),eye(6),inv(pluho([0.7822185856 0.5693963449 0.2528277810 -0.1671338270;-0.5693349884 0.8181046273 -0.0810091957 -0.0000077493;-0.2529659173 -0.0805768033 0.9641139058 0.0000038505;0.0000000000 0.0000000000 0.0000000000 1.0000000000])),eye(6),eye(6),inv(pluho([-0.2045694583 0.9759499257 0.0753198473 0.0000008108;-0.8631804239 -0.1435746457 -0.4840515231 0.0000069167;-0.4615960275 -0.1640367754 0.8717918580 -0.2885936396;0.0000000000 0.0000000000 0.0000000000 1.0000000000])),eye(6),inv(pluho([0.9282140708 0.2806752023 -0.2442131647 -0.0000090102;0.0276332465 0.6025814777 0.7975788153 0.0000013469;0.3710189249 -0.7470722813 0.5515686393 -0.2501673088;0.0000000000 0.0000000000 0.0000000000 1.0000000000])),eye(6),inv(pluho([0.9918872616 0.1263804939 0.0136978459 0.0777156893;-0.1236667821 0.9842667695 -0.1261960911 0.0258498624;-0.0294310589 0.1234783267 0.9919107397 -0.0818365185;0.0000000000 0.0000000000 0.0000000000 1.0000000000])),eye(6),eye(6),inv(pluho([-0.9758220275 0.1236145002 0.1802521172 -0.0000007438;0.1004456556 0.9860848778 -0.1324661614 -0.0000009629;-0.1941186252 -0.1111578561 -0.9746598843 -0.4166133153;0.0000000000 0.0000000000 0.0000000000 1.0000000000])),inv(pluho([-0.9505142897 -0.2969732696 0.0912658859 -0.0000731275;-0.1242159360 0.0940061044 -0.9877921105 0.0063187126;0.2847683024 -0.9502471937 -0.1262429583 0.4113371128;0.0000000000 0.0000000000 0.0000000000 1.0000000000])),eye(6),inv(pluho([0.9925597327 -0.1201365671 -0.0198086410 -0.0862986734;0.1173659308 0.9873091527 -0.1069853975 0.0336316042;0.0324101110 0.1038645380 0.9940632487 -0.1027320081;0.0000000000 0.0000000000 0.0000000000 1.0000000000])),eye(6),eye(6),inv(pluho([-0.9311665210 -0.1192092456 -0.3445548808 -0.0000020029;0.0179768967 0.9288754709 -0.3699556605 -0.0000007862;0.3641507122 -0.3506843528 -0.8627947284 -0.4197309109;0.0000000000 0.0000000000 0.0000000000 1.0000000000])),inv(pluho([-0.8553690247 0.5123695363 -0.0762973759 -0.0278996413;0.0050693626 -0.1390003416 -0.9902793578 0.1115025626;-0.5179943369 -0.8474410676 0.1162991997 0.3698874443;0.0000000000 0.0000000000 0.0000000000 1.0000000000])),eye(6)};
model.I = {mcI(1.00000,[0.0000000000 0.0000000000 0.0000000000],[1.0000000000 0.0000000000 0.0000000000;0.0000000000 1.0000000000 0.0000000000;0.0000000000 0.0000000000 1.0000000000]),zeros(6,6),zeros(6,6),mcI(1.00000,[0.0000000000 0.0000000000 0.0000000000],[1.0000000000 0.0000000000 0.0000000000;0.0000000000 1.0000000000 0.0000000000;0.0000000000 0.0000000000 1.0000000000]),zeros(6,6),zeros(6,6),mcI(1.00000,[0.0000000000 0.0000000000 0.0000000000],[1.0000000000 0.0000000000 0.0000000000;0.0000000000 1.0000000000 0.0000000000;0.0000000000 0.0000000000 1.0000000000]),zeros(6,6),mcI(1.00000,[0.0000000000 0.0000000000 0.0000000000],[1.0000000000 0.0000000000 0.0000000000;0.0000000000 1.0000000000 0.0000000000;0.0000000000 0.0000000000 1.0000000000]),zeros(6,6),zeros(6,6),mcI(1.00000,[0.0000000000 0.0000000000 0.0000000000],[1.0000000000 0.0000000000 0.0000000000;0.0000000000 1.0000000000 0.0000000000;0.0000000000 0.0000000000 1.0000000000]),zeros(6,6),mcI(1.00000,[0.0000000000 0.0000000000 0.0000000000],[1.0000000000 0.0000000000 0.0000000000;0.0000000000 1.0000000000 0.0000000000;0.0000000000 0.0000000000 1.0000000000]),zeros(6,6),mcI(1.00000,[0.0000000000 0.0000000000 0.0000000000],[1.0000000000 0.0000000000 0.0000000000;0.0000000000 1.0000000000 0.0000000000;0.0000000000 0.0000000000 1.0000000000]),zeros(6,6),mcI(1.00000,[0.0000000000 0.0000000000 0.0000000000],[1.0000000000 0.0000000000 0.0000000000;0.0000000000 1.0000000000 0.0000000000;0.0000000000 0.0000000000 1.0000000000]),zeros(6,6),zeros(6,6),mcI(1.00000,[0.0000000000 0.0000000000 0.0000000000],[1.0000000000 0.0000000000 0.0000000000;0.0000000000 1.0000000000 0.0000000000;0.0000000000 0.0000000000 1.0000000000]),zeros(6,6),mcI(1.00000,[0.0000000000 0.0000000000 0.0000000000],[1.0000000000 0.0000000000 0.0000000000;0.0000000000 1.0000000000 0.0000000000;0.0000000000 0.0000000000 1.0000000000]),zeros(6,6),mcI(1.00000,[0.0000000000 0.0000000000 0.0000000000],[1.0000000000 0.0000000000 0.0000000000;0.0000000000 1.0000000000 0.0000000000;0.0000000000 0.0000000000 1.0000000000]),zeros(6,6),zeros(6,6),mcI(1.00000,[0.0000000000 0.0000000000 0.0000000000],[1.0000000000 0.0000000000 0.0000000000;0.0000000000 1.0000000000 0.0000000000;0.0000000000 0.0000000000 1.0000000000]),mcI(1.00000,[0.0000000000 0.0000000000 0.0000000000],[1.0000000000 0.0000000000 0.0000000000;0.0000000000 1.0000000000 0.0000000000;0.0000000000 0.0000000000 1.0000000000]),zeros(6,6),mcI(1.00000,[0.0000000000 0.0000000000 0.0000000000],[1.0000000000 0.0000000000 0.0000000000;0.0000000000 1.0000000000 0.0000000000;0.0000000000 0.0000000000 1.0000000000]),zeros(6,6),zeros(6,6),mcI(1.00000,[0.0000000000 0.0000000000 0.0000000000],[1.0000000000 0.0000000000 0.0000000000;0.0000000000 1.0000000000 0.0000000000;0.0000000000 0.0000000000 1.0000000000]),mcI(1.00000,[0.0000000000 0.0000000000 0.0000000000],[1.0000000000 0.0000000000 0.0000000000;0.0000000000 1.0000000000 0.0000000000;0.0000000000 0.0000000000 1.0000000000]),zeros(6,6),mcI(1.00000,[0.0000000000 0.0000000000 0.0000000000],[1.0000000000 0.0000000000 0.0000000000;0.0000000000 1.0000000000 0.0000000000;0.0000000000 0.0000000000 1.0000000000])};

model.markers.name = {'EIASD','CID','EIPSD','EIPSG','CIG','EIASG','MANU','MIDSTERNUM','XIPHOIDE','C7','D3','D10','ZYGD','TEMPD','GLABELLE','TEMPG','ZYGG','CLAV1D','CLAV2D','CLAV3D','ACRANTD','ACRPOSTD','SCAPD','DELTD','BICEPSD','TRICEPSD','EPICOND','EPITROD','OLE1D','OLE2D','BRACHD','BRACHANTD','ABRAPOSTD','ABRASANTD','ULNAD','RADIUSD','METAC5D','METAC2D','MIDMETAC3D','CLAV1G','CLAV2G','CLAV3G','ACRANTG','ACRPOSTG','SCAPG','DELTG','BICEPSG','TRICEPSG','EPICONG','EPITROG','OLE1G','OLE2G','BRACHG','BRACHANTG','ABRAPOSTG','ABRANTG','ULNAG','RADIUSG','METAC5G','METAC2G','MIDMETAC3G','ISCHIO1D','TFLD','ISCHIO2D','CONDEXTD','CONDINTD','CRETED','JAMBLATD','TUBD','MALEXTD','MALINTD','CONDEXTD','CONDINTD','CALCD','MIDMETA4D','MIDMETA1D','SCAPHOIDED','METAT5D','METAT1D','ISCHIO1G','TFLG','ISCHIO2G','CONEXTG','CONDINTG','CRETEG','JAMBLATG','TUBG','ACHILLEG','MALEXTG','MALINTG','CONEXTG','CONDINTG','CALCG','MIDMETA4G','MIDMETA1G','SCAPHOIDEG','METAT5G','METAT1G'};

model.markers.parent = [6 6 6 6 6 6 9 9 9 9 9 9 12 12 12 12 12 14 14 14 14 14 14 17 17 17 17 17 19 19 19 19 19 19 19 19 21 21 21 23 23 23 23 23 23 26 26 26 26 26 28 28 28 28 28 28 28 28 30 30 30 33 33 33 33 33 34 34 34 34 34 34 34 36 36 36 36 36 36 39 39 39 39 39 40 40 40 40 40 40 40 40 42 42 42 42 42 42];
model.markers.coordinates = [0.1124034822 0.0912852857 0.0000000000;0.1552048181 -0.0014049505 0.0850283471;0.0453337603 -0.0943935354 -0.0040233244;-0.0453337603 -0.0937214520 0.0040233244;-0.1671288207 -0.0013192037 0.0862814291;-0.1124034822 0.0968297018 0.0000000000;-0.0024072283 0.0723418694 0.1807887663;0.0000506959 0.1088839241 0.0945860654;0.0020055425 0.1142918883 0.0375705328;0.0024072283 -0.0723418694 0.2412646498;0.0020259785 -0.1083956068 0.1046635181;0.0020055425 -0.1033939475 0.0283087269;0.0751090214 0.1128652059 -0.0112489317;0.0756662634 0.1341243745 0.0215404096;0.0006188041 0.1693068689 -0.0266566810;-0.0812132226 0.1216103472 0.0146833118;-0.0757653168 0.1128652059 -0.0112489317;-0.0139660587 0.0487953896 0.0174465962;0.0872118569 0.0191537966 0.0441104091;0.1316350339 0.0137057570 0.0563961693;0.1585149886 0.0233625396 0.0474761556;0.1849578402 -0.0438278243 0.0393156819;0.0697540625 -0.1128493743 0.0474761556;0.0358713371 0.0186603418 -0.1348721489;-0.0085077302 0.0600162995 -0.2360651210;0.0365404059 -0.0307655107 -0.1875100671;0.0344431267 -0.0100810985 -0.2796138147;-0.0507912682 -0.0100810985 -0.2802689229;-0.0334102694 0.0028952204 0.0060947628;-0.0378301347 -0.0024487897 -0.0406354937;-0.0164156137 -0.0353393407 -0.0657854056;0.0439824967 0.0403328037 -0.0693720368;0.0184460454 -0.0249684313 -0.2100719429;0.0124763910 0.0424376617 -0.2058738451;-0.0193711968 -0.0147510892 -0.2443806235;0.0293240790 -0.0147510892 -0.2474670977;-0.0142059560 -0.0363630862 -0.0568428284;0.0562332334 -0.0363630862 -0.0568428284;0.0255838836 -0.0346156093 -0.0344497047;0.0057687218 0.0360913543 0.0242442462;-0.0919090549 0.0097440845 0.0524033740;-0.1362587920 0.0118185471 0.0585311934;-0.1573748487 0.0141942224 0.0536728724;-0.1940054097 -0.0371574545 0.0411278861;-0.0791046241 -0.1189843486 0.0536728724;-0.0437672122 0.0009301612 -0.1155733049;-0.0242007843 0.0537708095 -0.2238933172;-0.0392939245 -0.0407838030 -0.1790241409;-0.0394242043 -0.0139490290 -0.2718897089;0.0418015420 -0.0139490290 -0.2788746107;0.0316557402 0.0074324600 0.0059810967;0.0394153590 -0.0041029227 -0.0429283706;0.0144145194 -0.0497588257 -0.0728566929;-0.0553399171 0.0146342006 -0.0510500089;-0.0084476682 -0.0321426233 -0.2147803382;-0.0200523797 0.0335626717 -0.2065219183;0.0230540677 -0.0236933033 -0.2377298533;-0.0207380412 -0.0236933033 -0.2466006491;0.0082544823 0.0538122877 -0.0465311225;-0.0629249053 0.0538122877 -0.0465311225;-0.0271839495 0.0277505879 -0.0414542613;0.0351370805 -0.0963775318 -0.2073997965;0.1003595132 -0.0228632286 -0.2616573348;0.0425431674 -0.0916732430 -0.3131338971;0.1016932662 -0.0016617420 -0.3715697017;-0.0118592924 -0.0016617420 -0.4135316555;-0.0356919295 0.0366366059 0.1778844357;-0.1104856711 -0.0095779406 0.1310933727;-0.0485809772 0.0405404590 0.0569218268;-0.0394537989 -0.0127825085 0.4100921549;0.0394537989 0.0152704108 0.3780471863;-0.1081453017 0.0059335018 -0.0253509264;0.0108072099 -0.0034455995 -0.0049229443;0.0040618234 -0.0404348499 0.0571969107;0.0431735889 0.0027218993 -0.0748477845;-0.0395950004 0.0046985762 -0.0785532289;-0.0046583716 0.0394922645 -0.0557282091;0.0464356866 -0.0000000000 -0.1178767987;-0.0464356866 -0.0000000000 -0.1318892736;-0.0215142431 -0.0977133804 -0.1912648136;-0.1099595019 -0.0116000447 -0.2380337400;-0.0428953838 -0.0856623352 -0.2933688315;-0.1140887787 0.0096794671 -0.3748253776;-0.0072219065 0.0096794671 -0.3936467109;0.0329731661 0.0882710942 0.1519855039;0.1022657175 0.0436033668 0.1509575426;0.0728022370 0.0545570360 0.0379772421;-0.0120226393 0.0151283775 0.2764040818;0.0120226393 0.1027856606 0.3725912764;-0.0660909538 0.1110715688 0.3380220854;0.1227635424 0.0068424693 -0.0030187077;0.0163992705 0.0006935469 -0.0235968815;-0.0149575062 -0.0279926025 0.0506783427;-0.0321128803 0.0070115789 -0.0917269419;0.0415025418 -0.0048437811 -0.0889878158;0.0046573553 0.0355937732 -0.0611707549;-0.0452432141 0.0000000000 -0.1246363777;0.0452432141 0.0000000000 -0.1391311320]';

model.appearance.base = { 'line', [1.1 0 0; 0 0 0; 0 1.1 0; 0 0 0; 0 0 1.1]};

model.appearance.body{1} = {'sphere', model.markers.coordinates(:, 1), 0.01,'sphere', model.markers.coordinates(:, 2), 0.01,'sphere', model.markers.coordinates(:, 3), 0.01,'sphere', model.markers.coordinates(:, 4), 0.01,'sphere', model.markers.coordinates(:, 5), 0.01,'sphere', model.markers.coordinates(:, 6), 0.01};
model.appearance.body{2} = {};
model.appearance.body{3} = {};
model.appearance.body{4} = {'sphere', model.markers.coordinates(:, 7), 0.01,'sphere', model.markers.coordinates(:, 8), 0.01,'sphere', model.markers.coordinates(:, 9), 0.01,'sphere', model.markers.coordinates(:, 10), 0.01,'sphere', model.markers.coordinates(:, 11), 0.01,'sphere', model.markers.coordinates(:, 12), 0.01};
model.appearance.body{5} = {};
model.appearance.body{6} = {};
model.appearance.body{7} = {'sphere', model.markers.coordinates(:, 13), 0.01,'sphere', model.markers.coordinates(:, 14), 0.01,'sphere', model.markers.coordinates(:, 15), 0.01,'sphere', model.markers.coordinates(:, 16), 0.01,'sphere', model.markers.coordinates(:, 17), 0.01};
model.appearance.body{8} = {};
model.appearance.body{9} = {'sphere', model.markers.coordinates(:, 18), 0.01,'sphere', model.markers.coordinates(:, 19), 0.01,'sphere', model.markers.coordinates(:, 20), 0.01,'sphere', model.markers.coordinates(:, 21), 0.01,'sphere', model.markers.coordinates(:, 22), 0.01,'sphere', model.markers.coordinates(:, 23), 0.01};
model.appearance.body{10} = {};
model.appearance.body{11} = {};
model.appearance.body{12} = {'sphere', model.markers.coordinates(:, 24), 0.01,'sphere', model.markers.coordinates(:, 25), 0.01,'sphere', model.markers.coordinates(:, 26), 0.01,'sphere', model.markers.coordinates(:, 27), 0.01,'sphere', model.markers.coordinates(:, 28), 0.01};
model.appearance.body{13} = {};
model.appearance.body{14} = {'sphere', model.markers.coordinates(:, 29), 0.01,'sphere', model.markers.coordinates(:, 30), 0.01,'sphere', model.markers.coordinates(:, 31), 0.01,'sphere', model.markers.coordinates(:, 32), 0.01,'sphere', model.markers.coordinates(:, 33), 0.01,'sphere', model.markers.coordinates(:, 34), 0.01,'sphere', model.markers.coordinates(:, 35), 0.01,'sphere', model.markers.coordinates(:, 36), 0.01};
model.appearance.body{15} = {};
model.appearance.body{16} = {'sphere', model.markers.coordinates(:, 37), 0.01,'sphere', model.markers.coordinates(:, 38), 0.01,'sphere', model.markers.coordinates(:, 39), 0.01};
model.appearance.body{17} = {};
model.appearance.body{18} = {'sphere', model.markers.coordinates(:, 40), 0.01,'sphere', model.markers.coordinates(:, 41), 0.01,'sphere', model.markers.coordinates(:, 42), 0.01,'sphere', model.markers.coordinates(:, 43), 0.01,'sphere', model.markers.coordinates(:, 44), 0.01,'sphere', model.markers.coordinates(:, 45), 0.01};
model.appearance.body{19} = {};
model.appearance.body{20} = {};
model.appearance.body{21} = {'sphere', model.markers.coordinates(:, 46), 0.01,'sphere', model.markers.coordinates(:, 47), 0.01,'sphere', model.markers.coordinates(:, 48), 0.01,'sphere', model.markers.coordinates(:, 49), 0.01,'sphere', model.markers.coordinates(:, 50), 0.01};
model.appearance.body{22} = {};
model.appearance.body{23} = {'sphere', model.markers.coordinates(:, 51), 0.01,'sphere', model.markers.coordinates(:, 52), 0.01,'sphere', model.markers.coordinates(:, 53), 0.01,'sphere', model.markers.coordinates(:, 54), 0.01,'sphere', model.markers.coordinates(:, 55), 0.01,'sphere', model.markers.coordinates(:, 56), 0.01,'sphere', model.markers.coordinates(:, 57), 0.01,'sphere', model.markers.coordinates(:, 58), 0.01};
model.appearance.body{24} = {};
model.appearance.body{25} = {'sphere', model.markers.coordinates(:, 59), 0.01,'sphere', model.markers.coordinates(:, 60), 0.01,'sphere', model.markers.coordinates(:, 61), 0.01};
model.appearance.body{26} = {};
model.appearance.body{27} = {};
model.appearance.body{28} = {'sphere', model.markers.coordinates(:, 62), 0.01,'sphere', model.markers.coordinates(:, 63), 0.01,'sphere', model.markers.coordinates(:, 64), 0.01,'sphere', model.markers.coordinates(:, 65), 0.01,'sphere', model.markers.coordinates(:, 66), 0.01};
model.appearance.body{29} = {'sphere', model.markers.coordinates(:, 67), 0.01,'sphere', model.markers.coordinates(:, 68), 0.01,'sphere', model.markers.coordinates(:, 69), 0.01,'sphere', model.markers.coordinates(:, 70), 0.01,'sphere', model.markers.coordinates(:, 71), 0.01,'sphere', model.markers.coordinates(:, 72), 0.01,'sphere', model.markers.coordinates(:, 73), 0.01};
model.appearance.body{30} = {};
model.appearance.body{31} = {'sphere', model.markers.coordinates(:, 74), 0.01,'sphere', model.markers.coordinates(:, 75), 0.01,'sphere', model.markers.coordinates(:, 76), 0.01,'sphere', model.markers.coordinates(:, 77), 0.01,'sphere', model.markers.coordinates(:, 78), 0.01,'sphere', model.markers.coordinates(:, 79), 0.01};
model.appearance.body{32} = {};
model.appearance.body{33} = {};
model.appearance.body{34} = {'sphere', model.markers.coordinates(:, 80), 0.01,'sphere', model.markers.coordinates(:, 81), 0.01,'sphere', model.markers.coordinates(:, 82), 0.01,'sphere', model.markers.coordinates(:, 83), 0.01,'sphere', model.markers.coordinates(:, 84), 0.01};
model.appearance.body{35} = {'sphere', model.markers.coordinates(:, 85), 0.01,'sphere', model.markers.coordinates(:, 86), 0.01,'sphere', model.markers.coordinates(:, 87), 0.01,'sphere', model.markers.coordinates(:, 88), 0.01,'sphere', model.markers.coordinates(:, 89), 0.01,'sphere', model.markers.coordinates(:, 90), 0.01,'sphere', model.markers.coordinates(:, 91), 0.01,'sphere', model.markers.coordinates(:, 92), 0.01};
model.appearance.body{36} = {};
model.appearance.body{37} = {'sphere', model.markers.coordinates(:, 93), 0.01,'sphere', model.markers.coordinates(:, 94), 0.01,'sphere', model.markers.coordinates(:, 95), 0.01,'sphere', model.markers.coordinates(:, 96), 0.01,'sphere', model.markers.coordinates(:, 97), 0.01,'sphere', model.markers.coordinates(:, 98), 0.01};
model = floatbase(model);

model.nq = model.NB;
model.nx = model.nq+model.nq;
model.nu = model.nq-6;

data.x0 = zeros(model.nx,1);
data.u0 = zeros(model.nu,1);

model.idx_q = 1:model.nq;
model.idx_v = model.nq+1:2*model.nq;

qmin_base = [-inf,-inf,-inf,-inf,-pi/4,-inf];
qmax_base = [ inf, inf, inf, inf, pi/4, inf];
qmin_thorax = [-pi/2,-pi/2.1,-pi/2];
qmax_thorax = [ pi/2, pi/2.1, pi/2];
qmin_tete = [-pi/2,-pi/2.1,-pi/2];
qmax_tete = [ pi/2, pi/2.1, pi/2];
qmin_epaule_droite = [-pi/2,-pi/2];
qmax_epaule_droite = [ pi/2, pi/2];
qmin_bras_droit = [-pi,-pi/2.1,-pi];
qmax_bras_droit = [ pi, pi/2.1, pi];
qmin_avantbras_droit = [ 0,-pi/2];
qmax_avantbras_droit = [pi, pi/2];
qmin_main_droite = [-pi/2,-pi/2];
qmax_main_droite = [ pi/2, pi/2];
qmin_epaule_gauche = [-pi/2,-pi/2];
qmax_epaule_gauche = [ pi/2, pi/2];
qmin_bras_gauche = [-pi,-pi/2.1,-pi];
qmax_bras_gauche = [ pi, pi/2.1, pi];
qmin_avantbras_gauche = [ 0,-pi/2];
qmax_avantbras_gauche = [pi, pi/2];
qmin_main_gauche = [-pi/2,-pi/2];
qmax_main_gauche = [ pi/2, pi/2];
qmin_cuisse_droite = [-pi,-pi/2.1,-pi/2];
qmax_cuisse_droite = [ pi, pi/2.1, pi/2];
qmin_jambe_droite = [-pi];
qmax_jambe_droite = [  0];
qmin_pied_droit = [-pi/2,-pi/2];
qmax_pied_droit = [ pi/2, pi/2];
qmin_cuisse_gauche = [-pi,-pi/2.1,-pi/2];
qmax_cuisse_gauche = [ pi, pi/2.1, pi/2];
qmin_jambe_gauche = [-pi];
qmax_jambe_gauche = [  0];
qmin_pied_gauche = [-pi/2,-pi/2];
qmax_pied_gauche = [ pi/2, pi/2];

qdotmin_base = [-inf,-inf,-inf,-inf,-inf,-inf];
qdotmax_base = [ inf, inf, inf, inf, inf, inf];

model.xmin = [qmin_base'; ... % q
	qmin_thorax'; ...
	qmin_tete'; ...
	qmin_epaule_droite'; ...
	qmin_bras_droit'; ...
	qmin_avantbras_droit'; ...
	qmin_main_droite'; ...
	qmin_epaule_gauche'; ...
	qmin_bras_gauche'; ...
	qmin_avantbras_gauche'; ...
	qmin_main_gauche'; ...
	qmin_cuisse_droite'; ...
	qmin_jambe_droite'; ...
	qmin_pied_droit'; ...
	qmin_cuisse_gauche'; ...
	qmin_jambe_gauche'; ...
	qmin_pied_gauche'; ...
	qdotmin_base'; -200*ones(model.nq-6,1)]; % qdot
model.xmax = [qmax_base';  ... % q
	qmax_thorax'; ...
	qmax_tete'; ...
	qmax_epaule_droite'; ...
	qmax_bras_droit'; ...
	qmax_avantbras_droit'; ...
	qmax_main_droite'; ...
	qmax_epaule_gauche'; ...
	qmax_bras_gauche'; ...
	qmax_avantbras_gauche'; ...
	qmax_main_gauche'; ...
	qmax_cuisse_droite'; ...
	qmax_jambe_droite'; ...
	qmax_pied_droit'; ...
	qmax_cuisse_gauche'; ...
	qmax_jambe_gauche'; ...
	qmax_pied_gauche'; ...
	qdotmax_base';  200*ones(model.nq-6,1)]; % qdot

model.umin = -150*ones(model.nu,1);
model.umax =  150*ones(model.nu,1);

if isfield(data, 'gravity')
	model.gravity = data.gravity;
end
end
