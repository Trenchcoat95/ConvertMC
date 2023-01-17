    #!/bin/bash
    root.exe -n -b -l <<\EOF 2>&1 | tee makeData.log
    gSystem->AddIncludePath("-I\"/home/federico/Documents/Universita/Federico_2020-2021/Aliwork/fastMCKalman/fastMCKalman/aliKalman/test/\"")
    gSystem->AddIncludePath("-I\"/home/federico/Documents/Universita/Federico_2020-2021/Aliwork/fastMCKalman/fastMCKalman/MC/\"")
    gSystem->Load("/home/federico/Documents/Universita/Federico_2020-2021/Aliwork/fastMCKalman/fastMCKalman/aliKalman/test/AliExternalTrackParam.so");
    .L /home/federico/Documents/Universita/Federico_2020-2021/Aliwork/fastMCKalman/fastMCKalman/MC/fastSimulation.cxx++g
    .L /home/federico/Documents/Universita/Federico_2020-2021/Aliwork/ConvertMC/testNDGAr.C++g
    AliPDG::AddParticlesToPdgDataBase();
    AliLog::SetPrintRepetitions(0);
    testNDGAr(10,kTRUE);            //setup for the looper development
    .q
