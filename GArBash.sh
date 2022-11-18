
# source $fastMCKalman/fastMCKalman/tests/unitTest.sh


alias helpCat=cat
[[ -x "$(command -v pygmentize)" ]] && alias helpCat="pygmentize -O style=borland,linenos=1 -l bash"
init(){
  cat <<HELP_USAGE | helpCat
  makeData
  analyzeLogs
HELP_USAGE
}

makeData(){
     export nPoints=${1:-40000}

    cat <<EOF >  makeData.sh
    #!/bin/bash
    root.exe -n -b -l <<\EOF 2>&1 | tee makeData.log
    gSystem->AddIncludePath("-I\"$fastMCKalman/fastMCKalman/aliKalman/test/\"")
    gSystem->AddIncludePath("-I\"$fastMCKalman/fastMCKalman/MC/\"")
    gSystem->Load("$fastMCKalman/fastMCKalman/aliKalman/test/AliExternalTrackParam.so");
    .L $fastMCKalman/fastMCKalman/MC/fastSimulation.cxx++g
    .L $ConvertMC/testNDGAr.C++g
    AliPDG::AddParticlesToPdgDataBase();
    AliLog::SetPrintRepetitions(0);
    testNDGAr(${nPoints},kTRUE);            //setup for the looper development
    .q
EOF
   chmod a+x makeData.sh

   ./makeData.sh
}


analyzeLogs(){
   errors=( "short track" "Too few consecutive points" "Rotation failed" "Propagation failed" "Update failed" "Too big chi2"
   "Correct for material failed" "PropagateToMirrorX failed" )
   errorSources=( "fastParticle::reconstructParticleFull:")
   for errorSource in  "${errorSources[@]}"; do
      nErrorstot=$(cat makeData.log | grep -c ${errorSource})
      echo  ${errorSource} ${nErrorstot}
      for error in "${errors[@]}"; do
        nErrors=$(cat makeData.log | grep ${errorSource} | grep -c "${error}")
        echo ${errorSource} ${error} ${nErrors}
      done;
   done
}


init