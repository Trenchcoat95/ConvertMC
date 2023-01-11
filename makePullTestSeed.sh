#!/bin/bash
    echo fastParticle.root >fastParticle.list
    root.exe -n -b -l <<\EOF 2>&1 | tee makePullTestSeed.log
    gSystem->AddIncludePath("-I\"/home/federico/Documents/Universita/Federico_2020-2021/Aliwork/fastMCKalman/fastMCKalman/fastMCKalman/aliKalman/test/\"")
    gSystem->Load("/home/federico/Documents/Universita/Federico_2020-2021/Aliwork/fastMCKalman/fastMCKalman/aliKalman/test/AliExternalTrackParam.so");
    .L /home/federico/Documents/Universita/Federico_2020-2021/Aliwork/fastMCKalman/fastMCKalman/MC/fastSimulation.cxx+g
    .L /home/federico/Documents/Universita/Federico_2020-2021/Aliwork/fastMCKalman/fastMCKalman/MC/fastSimulationTest.C+g
     .L /home/federico/Documents/Universita/Federico_2020-2021/Aliwork/fastMCKalman/fastMCKalman/MC/test_fastSimulation.C+g
     initTreeFast()
     testPullsSeed()
    .q
