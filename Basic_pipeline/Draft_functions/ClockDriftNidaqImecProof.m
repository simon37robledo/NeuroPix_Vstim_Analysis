squareWNi = load('PV103_Experiment_12_6_23_6_g0_tcat.nidq.xd_16_0_500.txt');
squareWap = load('PV103_Experiment_12_6_23_6_g0_tcat.imec0.ap.xd_384_6_500.txt');

figure;
plot(squareWNi,squareWap);hold on; plot(squareWNi,squareWNi)

squareWNi(end)

squareWap(end)



exnoS = load('PV103_Experiment_12_6_23_6_g0_tcat.nidq.xid_16_7_0.txt');

exS = load('out_7inv.txt');

(exnoS(end)-exS(end))*1000-(exnoS(1)-exS(1))*1000

(exnoS(1)-exS(1))*1000