function [model, data] = GenerateModel(data)
% version 1

% File extracted from DoCi.s2mMod

% Informations générales
% root_actuated	1
% external_forces	0

import casadi.*

model.name = 'DoCi';

model.NB = 37;

model.bodyN_name = {'Tx_pelvis','Ty_pelvis','Tz_pelvis','Rx_pelvis','Ry_pelvis','Rz_pelvis','Rx_thorax','Ry_thorax','Rz_thorax','Rx_tete','Ry_tete','Rz_tete','Ry_epauled','Rz_epauled','Rx_brasd','Ry_brasd','Rz_brasd','Rx_abrasd','Rz_abrasd','Rx_maind','Ry_maind','Ry_epauleg','Rz_epauleg','Rx_brasg','Ry_brasg','Rz_brasg','Rx_abrasg','Rz_abrasg','Rx_maing','Ry_maing','Rx_cuissed','Ry_cuissed','Rz_cuissed','Rx_jambed','Rx_piedd','Rz_piedd','Rx_cuisseg','Ry_cuisseg','Rz_cuisseg','Rx_jambeg','Rx_piedg','Rz_piedg'};
model.jtype = {'R','Rx','Ry','Rz','Rx','Ry','Rz','Ry','Rz','Rx','Ry','Rz','Rx','Rz','Rx','Ry','Ry','Rz','Rx','Ry','Rz','Rx','Rz','Rx','Ry','Rx','Ry','Rz','Rx','Rx','Rz','Rx','Ry','Rz','Rx','Rx','Rz'};
model.parent = [0,1,2,3,4,5,6,4,8,9,10,11,12,13,14,15,4,17,18,19,20,21,22,23,24,1,26,27,28,29,30,1,32,33,34,35,36];

model.Xtree = {eye(6),inv(pluho([0.9978916035 0.0226102511 0.0608368653 -0.0018212712;-0.0082063648 0.9737889385 -0.2273054305 0.0138886119;-0.0643816993 0.2263269310 0.9719213533 0.1596893073;0.0000000000 0.0000000000 0.0000000000 1.0000000000])),eye(6),eye(6),inv(pluho([1 0 0 -0.0091056578;0 1 0 0.0110294186;0 0 1 0.3211734003;0 0 0 1])),eye(6),eye(6),inv(pluho([0.9730592445 0.2081178319 0.0992102548 0.0400367266;-0.2105715447 0.9774660818 0.0148217192 0.0168400533;-0.0938899949 -0.0353132675 0.9949561004 0.2064402092;0.0000000000 0.0000000000 0.0000000000 1.0000000000])),eye(6),inv(pluho([0.8520926065 0.3372449655 -0.4002549474 0.1336035708;-0.3959352037 0.9154883932 -0.0715284303 -0.0000016592;0.3423061559 0.2194238709 0.9136080452 0.0000052680;0.0000000000 0.0000000000 0.0000000000 1.0000000000])),eye(6),eye(6),inv(pluho([0.7283551265 -0.5976584780 -0.3351166266 0.0000177437;0.5344502598 0.8015893914 -0.2679876255 0.0000127460;0.4287910092 0.0160869929 0.9032604712 -0.2572510774;0.0000000000 0.0000000000 0.0000000000 1.0000000000])),eye(6),inv(pluho([0.9695766286 0.2442431497 -0.0163231515 0.0000002970;-0.0202423239 0.0135445967 -0.9997033522 0.0000045305;-0.2439496050 0.9696194241 0.0180765731 -0.2208340024;0.0000000000 0.0000000000 0.0000000000 1.0000000000])),eye(6),inv(pluho([0.9858361205 -0.1657087506 -0.0258409271 -0.0430360037;0.1616984529 0.9800283168 -0.1157501986 0.0117997144;0.0445056611 0.1099322888 0.9929421624 0.1973975518;0.0000000000 0.0000000000 0.0000000000 1.0000000000])),eye(6),inv(pluho([0.8990763455 -0.3633699329 0.2441802959 -0.1377749097;0.3567609066 0.9313849447 0.0724136752 -0.0000065312;-0.2537388037 0.0220085613 0.9670223590 -0.0000009607;0.0000000000 0.0000000000 0.0000000000 1.0000000000])),eye(6),eye(6),inv(pluho([0.7242952314 0.5604424486 0.4016225586 -0.0000114269;-0.4726994156 0.8276765763 -0.3025008227 0.0000063257;-0.5019478860 0.0292531546 0.8644030153 -0.2561324864;0.0000000000 0.0000000000 0.0000000000 1.0000000000])),eye(6),inv(pluho([0.8815292732 -0.4226194185 -0.2104731991 -0.0000070224;0.0202751852 0.4792710284 -0.8774327316 0.0000140342;0.4716938174 0.7692152549 0.4310601284 -0.2192288097;0.0000000000 0.0000000000 0.0000000000 1.0000000000])),eye(6),inv(pluho([0.9977104468 -0.0095289495 -0.0669556827 0.0815784902;0.0093853078 0.9999529322 -0.0024595573 0.0241591400;0.0669759682 0.0018255263 0.9977529189 -0.1080733352;0.0000000000 0.0000000000 0.0000000000 1.0000000000])),eye(6),eye(6),inv(pluho([0.9921459658 0.0061351118 -0.1249349552 0.0000119089;-0.0191894986 0.9944384703 -0.1035562255 -0.0000039273;0.1236047967 0.1051403305 0.9867458463 -0.3633578278;0.0000000000 0.0000000000 0.0000000000 1.0000000000])),inv(pluho([0.9392664688 0.3134751551 -0.1396847448 -0.0312824628;-0.1813709936 0.1078766364 -0.9774800222 -0.0000015767;-0.2913469813 0.9434489698 0.1581802076 -0.3512915783;0.0000000000 0.0000000000 0.0000000000 1.0000000000])),eye(6),inv(pluho([0.9935221802 0.0095973183 0.1132323674 -0.1009874211;-0.0141819469 0.9991088934 0.0397528800 0.0290780961;-0.1127499443 -0.0411012234 0.9927729547 -0.1018191597;0.0000000000 0.0000000000 0.0000000000 1.0000000000])),eye(6),eye(6),inv(pluho([0.9817948019 -0.0666617644 0.1778628020 0.0000041108;0.1063199717 0.9688250179 -0.2237725368 0.0000016957;-0.1574008603 0.2386090815 0.9582748435 -0.3651353441;0.0000000000 0.0000000000 0.0000000000 1.0000000000])),inv(pluho([0.9174311552 -0.3726105171 0.1395760661 0.0308964960;0.2434286916 0.2481321691 -0.9376422018 0.0001551783;0.3147420336 0.8941989876 0.3183482732 -0.3415499844;0.0000000000 0.0000000000 0.0000000000 1.0000000000])),eye(6)};
model.I = {mcI(7.03000,[0.0000000000 0.0000000000 0.0656900000],[0.0268500000 0.0000000000 0.0000000000;0.0000000000 0.0517500000 0.0000000000;0.0000000000 0.0000000000 0.0618400000]),zeros(6,6),zeros(6,6),mcI(16.34000,[0.0000000000 0.0000000000 0.1721000000],[0.3609000000 0.0000000000 0.0000000000;0.0000000000 0.4135000000 0.0000000000;0.0000000000 0.0000000000 0.1664000000]),zeros(6,6),zeros(6,6),mcI(4.18000,[0.0000000000 0.0000000000 0.0734000000],[0.0554000000 0.0000000000 0.0000000000;0.0000000000 0.0554000000 0.0000000000;0.0000000000 0.0000000000 0.0119000000]),zeros(6,6),mcI(0.78000,[0.1032000000 0.0000000000 0.0000000000],[0.0000000000 0.0000000000 0.0000000000;0.0000000000 0.0000000000 0.0000000000;0.0000000000 0.0000000000 0.0000000000]),zeros(6,6),zeros(6,6),mcI(1.54000,[0.0000000000 0.0000000000 -0.1052100000],[0.0069200000 0.0000000000 0.0000000000;0.0000000000 0.0069200000 0.0000000000;0.0000000000 0.0000000000 0.0017600000]),zeros(6,6),mcI(0.84000,[0.0000000000 0.0000000000 -0.1051000000],[0.0036000000 0.0000000000 0.0000000000;0.0000000000 0.0036000000 0.0000000000;0.0000000000 0.0000000000 0.0005000000]),zeros(6,6),mcI(0.33000,[0.0193500401 -0.0386228178 0.0303525541],[0.0010000000 0.0000000000 0.0000000000;0.0000000000 0.0011000000 0.0000000000;0.0000000000 0.0000000000 0.0001000000]),zeros(6,6),mcI(0.78000,[-0.1032000000 0.0000000000 0.0000000000],[0.0000000000 0.0000000000 0.0000000000;0.0000000000 0.0000000000 0.0000000000;0.0000000000 0.0000000000 0.0000000000]),zeros(6,6),zeros(6,6),mcI(1.54000,[0.0000000000 0.0000000000 -0.1052100000],[0.0069200000 0.0000000000 0.0000000000;0.0000000000 0.0069200000 0.0000000000;0.0000000000 0.0000000000 0.0017600000]),zeros(6,6),mcI(0.84000,[0.0000000000 0.0000000000 -0.1051000000],[0.0036000000 0.0000000000 0.0000000000;0.0000000000 0.0036000000 0.0000000000;0.0000000000 0.0000000000 0.0005000000]),zeros(6,6),mcI(0.33000,[-0.0268379198 -0.0479202741 0.0160504777],[0.0010000000 0.0000000000 0.0000000000;0.0000000000 0.0011000000 0.0000000000;0.0000000000 0.0000000000 0.0001000000]),zeros(6,6),zeros(6,6),mcI(8.71000,[0.0000000000 0.0000000000 -0.1646700000],[0.1136300000 0.0000000000 0.0000000000;0.0000000000 0.1136300000 0.0000000000;0.0000000000 0.0000000000 0.0351400000]),mcI(3.22000,[0.0000000000 0.0000000000 -0.1546000000],[0.0391000000 0.0000000000 0.0000000000;0.0000000000 0.0391000000 0.0000000000;0.0000000000 0.0000000000 0.0047000000]),zeros(6,6),mcI(0.83000,[0.0000000000 0.0000000000 -0.0719000000],[0.0052000000 0.0000000000 0.0000000000;0.0000000000 0.0051000000 0.0000000000;0.0000000000 0.0000000000 0.0007000000]),zeros(6,6),zeros(6,6),mcI(8.71000,[0.0000000000 0.0000000000 -0.1646700000],[0.1136300000 0.0000000000 0.0000000000;0.0000000000 0.1136300000 0.0000000000;0.0000000000 0.0000000000 0.0351400000]),mcI(3.22000,[0.0000000000 0.0000000000 -0.1546000000],[0.0391000000 0.0000000000 0.0000000000;0.0000000000 0.0391000000 0.0000000000;0.0000000000 0.0000000000 0.0047000000]),zeros(6,6),mcI(0.83000,[0.0000000000 0.0000000000 -0.0719000000],[0.0052000000 0.0000000000 0.0000000000;0.0000000000 0.0051000000 0.0000000000;0.0000000000 0.0000000000 0.0007000000])};

model.markers.name = {'EIASD','CID','EIPSD','EIPSG','CIG','EIASG','MANU','MIDSTERNUM','XIPHOIDE','C7','D3','D10','ZYGD','TEMPD','GLABELLE','TEMPG','ZYGG','CLAV1D','CLAV2D','CLAV3D','ACRANTD','ACRPOSTD','SCAPD','DELTD','BICEPSD','TRICEPSD','EPICOND','EPITROD','OLE1D','OLE2D','BRACHD','BRACHANTD','ABRAPOSTD','ABRASANTD','ULNAD','RADIUSD','METAC5D','METAC2D','MIDMETAC3D','CLAV1G','CLAV2G','CLAV3G','ACRANTG','ACRPOSTG','SCAPG','DELTG','BICEPSG','TRICEPSG','EPICONG','EPITROG','OLE1G','OLE2G','BRACHG','BRACHANTG','ABRAPOSTG','ABRANTG','ULNAG','RADIUSG','METAC5G','METAC2G','MIDMETAC3G','ISCHIO1D','TFLD','ISCHIO2D','CONDEXTD','CONDINTD','CRETED','JAMBLATD','TUBD','ACHILED','MALEXTD','MALINTD','CALCD','MIDMETA4D','MIDMETA1D','SCAPHOIDED','METAT5D','METAT1D','ISCHIO1G','TFLG','ISCHIO2G','CONEXTG','CONDINTG','CRETEG','JAMBLATG','TUBG','ACHILLEG','MALEXTG','MALINTG','CALCG','MIDMETA4G','MIDMETA1G','SCAPHOIDEG','METAT5G','METAT1G'};

model.markers.parent = [6 6 6 6 6 6 9 9 9 9 9 9 12 12 12 12 12 14 14 14 14 14 14 17 17 17 17 17 19 19 19 19 19 19 19 19 21 21 21 23 23 23 23 23 23 26 26 26 26 26 28 28 28 28 28 28 28 28 30 30 30 33 33 33 33 33 34 34 34 34 34 34 36 36 36 36 36 36 39 39 39 39 39 40 40 40 40 40 40 42 42 42 42 42 42];
model.markers.coordinates = [0.1281892217 0.0749765205 -0.0000000000;0.1283565815 0.0598085539 0.0840360207;0.0463989796 -0.0870102025 -0.0020821297;-0.0463989796 -0.0818536384 0.0020821297;-0.1331230125 0.0624886367 0.0986735213;-0.1281892217 0.0938873204 -0.0000000000;0.0033092731 0.0583861677 0.2185593169;-0.0001487671 0.0962082313 0.1543837235;0.0023228674 0.1049518289 0.0063884995;-0.0033092731 -0.0583861677 0.2887522678;-0.0006155205 -0.1020477242 0.1483691062;0.0023228674 -0.1122597960 -0.0029449590;0.0543 0.0964 0.0525;0.0544 0.0931 0.0845;-0.0181 0.1343 0.0708;-0.0815 0.0798 0.0684;-0.0741 0.0870 0.0363;-0.0204251174 0.0347019537 0.0251325052;0.0395918602 0.0331820039 0.0271597446;0.0900578780 0.0148499413 0.0525999129;0.1243493578 0.0237700044 0.0486119089;0.1443451425 -0.0358938126 0.0491604305;0.0645887672 -0.0966185426 0.0486119089;0.0294683691 0.0232712785 -0.1160801414;-0.0031585598 0.0509295887 -0.1872410486;0.0466575913 -0.0251498276 -0.1585434582;0.0408710889 -0.0115870383 -0.2417610291;-0.0393522180 -0.0115870383 -0.2413461931;-0.0233422978 -0.0211055800 0.0050432603;-0.0207958031 -0.0262882075 -0.0418521826;0.0262970393 -0.0384124900 -0.0524310731;0.0320069480 0.0405839563 -0.0697367917;0.0141429993 -0.0295895718 -0.1765179885;0.0068642577 0.0364858203 -0.1713768662;-0.0201238545 -0.0209682933 -0.2061854256;0.0283177194 -0.0209682933 -0.2100530340;-0.0099786710 -0.0546490682 0.0252250353;0.0546575097 -0.0546490682 0.0252250353;0.0193500401 -0.0386228178 0.0303525541;0.0214366534 0.0433643166 0.0310098816;-0.0470410398 0.0370749835 0.0290441687;-0.0909742198 0.0194200019 0.0492879390;-0.1223973823 0.0255765273 0.0446424337;-0.1463756201 -0.0371625636 0.0424684629;-0.0588404448 -0.0963499824 0.0446424337;-0.0374801686 0.0103542301 -0.1100486773;-0.0220222336 0.0455813520 -0.1749210879;-0.0463370944 -0.0279672066 -0.1641867303;-0.0324361317 -0.0166352200 -0.2423834441;0.0354166626 -0.0166352200 -0.2489670325;0.0249508783 -0.0189659871 0.0025586539;0.0218152196 -0.0262167646 -0.0453006319;-0.0273982356 -0.0394080794 -0.0659806015;-0.0291230058 0.0429020240 -0.0532180856;-0.0160561475 -0.0340425600 -0.1713610255;0.0016600856 0.0349155973 -0.1714542778;0.0166358836 -0.0264409289 -0.2012023595;-0.0238884713 -0.0264409289 -0.2067150586;0.0032823015 -0.0663103784 0.0047278884;-0.0595883318 -0.0663103784 0.0047278884;-0.0268379198 -0.0479202741 0.0160504777;-0.0161222809 -0.0868954611 -0.1438549586;0.0858782553 -0.0355421041 -0.1879816336;-0.0052912137 -0.0839930995 -0.2225152234;0.0803215448 0.0131698948 -0.3302596721;-0.0438125857 0.0131698948 -0.3412502486;-0.0055344456 0.0343440339 -0.1582887419;0.0809293631 -0.0092206446 -0.1276891350;0.0320599183 0.0484826103 -0.0473513964;-0.0109846796 -0.0608102118 -0.2774813224;0.0109846796 -0.0004928281 -0.3372788913;-0.0697005520 0.0026862504 -0.3234536431;0.0254033951 -0.0280072063 0.0495738746;0.0248382816 0.0087328731 -0.0753814697;-0.0437841517 0.0015823130 -0.0755595971;-0.0135895007 0.0307625535 -0.0589126193;0.0469196283 0.0000000000 -0.1059600338;-0.0469196283 0.0000000000 -0.1220510997;0.0032495654 -0.0873090508 -0.1395860133;-0.0885049688 -0.0220204867 -0.1672545657;0.0021047949 -0.0834387970 -0.2150065897;-0.0740478780 0.0120646217 -0.3238347560;0.0521425026 0.0120646217 -0.3326410093;0.0095817707 0.0300399852 -0.1611679064;-0.0697915784 -0.0059357757 -0.1346084865;-0.0110461889 0.0439587874 -0.0560266998;0.0061164913 -0.0705022425 -0.2529768558;-0.0061164913 -0.0132728309 -0.3369751028;0.0729537684 -0.0052393831 -0.3181682493;-0.0287623600 -0.0401323493 0.0570282373;-0.0252494126 0.0061766839 -0.0713369752;0.0443572046 -0.0065159690 -0.0702614227;0.0093760024 0.0266594901 -0.0463849087;-0.0465342877 0.0000000000 -0.1003520007;0.0465342877 0.0000000000 -0.1126042619]';

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
model.appearance.body{29} = {'sphere', model.markers.coordinates(:, 67), 0.01,'sphere', model.markers.coordinates(:, 68), 0.01,'sphere', model.markers.coordinates(:, 69), 0.01,'sphere', model.markers.coordinates(:, 70), 0.01,'sphere', model.markers.coordinates(:, 71), 0.01,'sphere', model.markers.coordinates(:, 72), 0.01};
model.appearance.body{30} = {};
model.appearance.body{31} = {'sphere', model.markers.coordinates(:, 73), 0.01,'sphere', model.markers.coordinates(:, 74), 0.01,'sphere', model.markers.coordinates(:, 75), 0.01,'sphere', model.markers.coordinates(:, 76), 0.01,'sphere', model.markers.coordinates(:, 77), 0.01,'sphere', model.markers.coordinates(:, 78), 0.01};
model.appearance.body{32} = {};
model.appearance.body{33} = {};
model.appearance.body{34} = {'sphere', model.markers.coordinates(:, 79), 0.01,'sphere', model.markers.coordinates(:, 80), 0.01,'sphere', model.markers.coordinates(:, 81), 0.01,'sphere', model.markers.coordinates(:, 82), 0.01,'sphere', model.markers.coordinates(:, 83), 0.01};
model.appearance.body{35} = {'sphere', model.markers.coordinates(:, 84), 0.01,'sphere', model.markers.coordinates(:, 85), 0.01,'sphere', model.markers.coordinates(:, 86), 0.01,'sphere', model.markers.coordinates(:, 87), 0.01,'sphere', model.markers.coordinates(:, 88), 0.01,'sphere', model.markers.coordinates(:, 89), 0.01};
model.appearance.body{36} = {};
model.appearance.body{37} = {'sphere', model.markers.coordinates(:, 90), 0.01,'sphere', model.markers.coordinates(:, 91), 0.01,'sphere', model.markers.coordinates(:, 92), 0.01,'sphere', model.markers.coordinates(:, 93), 0.01,'sphere', model.markers.coordinates(:, 94), 0.01,'sphere', model.markers.coordinates(:, 95), 0.01};
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
qmin_thorax = [-pi/2,-pi/2,-pi/2];
qmax_thorax = [ pi/2, pi/2, pi/2];
qmin_tete = [-pi/2,-pi/2,-pi/2];
qmax_tete = [ pi/2, pi/2, pi/2];
qmin_epaule_droite = [-pi/2,-pi/2];
qmax_epaule_droite = [ pi/2, pi/2];
qmin_bras_droit = [-pi,-pi/2,-pi];
qmax_bras_droit = [ pi, pi/2, pi];
qmin_avantbras_droit = [ 0,-pi/2];
qmax_avantbras_droit = [pi, pi/2];
qmin_main_droite = [-pi/2,-pi/2];
qmax_main_droite = [ pi/2, pi/2];
qmin_epaule_gauche = [-pi/2,-pi/2];
qmax_epaule_gauche = [ pi/2, pi/2];
qmin_bras_gauche = [-pi,-pi/2,-pi];
qmax_bras_gauche = [ pi, pi/2, pi];
qmin_avantbras_gauche = [ 0,-pi/2];
qmax_avantbras_gauche = [pi, pi/2];
qmin_main_gauche = [-pi/2,-pi/2];
qmax_main_gauche = [ pi/2, pi/2];
qmin_cuisse_droite = [-pi,-pi/2,-pi/2];
qmax_cuisse_droite = [ pi, pi/2, pi/2];
qmin_jambe_droite = [-pi];
qmax_jambe_droite = [  0];
qmin_pied_droit = [-pi/2,-pi/2];
qmax_pied_droit = [ pi/2, pi/2];
qmin_cuisse_gauche = [-pi,-pi/2,-pi/2];
qmax_cuisse_gauche = [ pi, pi/2, pi/2];
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
              qdotmin_base'; -100*ones(model.nq-6,1)]; % qdot
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
              qdotmax_base';  100*ones(model.nq-6,1)]; % qdot

model.umin = -50*ones(model.nu,1);
model.umax =  50*ones(model.nu,1);

% model.gravity = [0 0 0];
end
