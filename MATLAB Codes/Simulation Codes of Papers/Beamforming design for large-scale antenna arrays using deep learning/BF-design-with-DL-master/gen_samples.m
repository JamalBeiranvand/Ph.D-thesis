
    pnr_dB = 20;
    ITER =10000;
    Lest =3;
    Num_paths =3;
    %gen pcsi and ecsi. 
    %pcsi: perfect csi ecsi:estimated csi
    [pcsi,  ecsi] = channel_gen_LOS(pnr_dB, ITER, Num_paths, Lest,1,64,1);
    save('pcsi.mat', 'pcsi')
    save('ecsi.mat', 'ecsi')
