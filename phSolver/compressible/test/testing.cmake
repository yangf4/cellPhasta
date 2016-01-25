macro(c_parallel_test name procs dir exe)
  set(tname compressible_${name})
  add_test(
    NAME ${tname}
    COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} ${procs} ${exe} ${ARGN}
    WORKING_DIRECTORY ${dir} )
  set_tests_properties(${tname} PROPERTIES LABELS "phsolver_compressible")
endmacro(c_parallel_test)

macro(c_serial_test name dir exe)
  set(tname compressible_${name})
  add_test( NAME ${tname} COMMAND ${exe} ${ARGN} WORKING_DIRECTORY ${dir})
  set_tests_properties(${tname} PROPERTIES LABELS "phsolver_compressible")
endmacro(c_serial_test)

set(CDIR ${CASES}/compressible_dg_interface)
c_serial_test(inpCfg_dg_interface ${CDIR} cp ${PHASTA_SOURCE_DIR}/phSolver/common/input.config ${CDIR})

c_serial_test(resetNumStart_dg_interface ${CDIR}
  cp ${CDIR}/numstart.dat ${CDIR}/1-procs_case/numstart.dat)
c_serial_test( dg_interface ${CDIR}/1-procs_case ${PHASTA_BINARY_DIR}/bin/phastaC.exe)

set(CDIR ${CASES}/compressible)
c_serial_test(inpCfg ${CDIR} cp ${PHASTA_SOURCE_DIR}/phSolver/common/input.config ${CDIR})

c_serial_test(linkProcsDir-sync ${CDIR}
  ln -snf ${CDIR}/2-procs_case-SyncIO-1 ${CDIR}/2-procs_case)
if(HAS_VALGRIND)
  c_serial_test(resetNumStartValgrind-sync ${CDIR}
    cp ${CDIR}/numstart.dat ${CDIR}/2-procs_case/numstart.dat)
  c_parallel_test(valgrind-sync 2 ${CDIR}
    valgrind --leak-check=yes --log-file=cSyncValgrind.%p
    ${PHASTA_BINARY_DIR}/bin/phastaC.exe)
endif(HAS_VALGRIND)
c_serial_test(resetNumStart-sync ${CDIR}
  cp ${CDIR}/numstart.dat ${CDIR}/2-procs_case/numstart.dat)
c_parallel_test( sync 2 ${CDIR} ${PHASTA_BINARY_DIR}/bin/phastaC.exe)
set(compareArgs
  ${CDIR}/2-procs_case-SyncIO-1/
  ${CDIR}/2-procs_case-SyncIO-1_ref/
  1 1e-6)
c_parallel_test(compare-sync 2 ${CDIR}
  ${PHASTA_BINARY_DIR}/bin/checkphasta ${compareArgs})
c_parallel_test(restart-sync 2 ${CDIR} ${PHASTA_BINARY_DIR}/bin/phastaC.exe)
c_parallel_test(compareRestart-sync 2 ${CDIR}
  ${PHASTA_BINARY_DIR}/bin/checkphasta ${compareArgs})

c_serial_test(linkProcsDir-posix ${CDIR}
  ln -snf ${CDIR}/2-procs_case-Posix ${CDIR}/2-procs_case)
if(HAS_VALGRIND)
  c_serial_test(resetNumStartValgrind-posix ${CDIR}
    cp ${CDIR}/numstart.dat ${CDIR}/2-procs_case/numstart.dat)
  c_parallel_test(compressibleValgrind-posix 2 ${CDIR}
    valgrind --leak-check=yes --log-file=cPosixValgrind.%p
    ${PHASTA_BINARY_DIR}/bin/phastaC.exe)
endif(HAS_VALGRIND)
c_serial_test(resetNumStart-posix ${CDIR}
  cp ${CDIR}/numstart.dat ${CDIR}/2-procs_case/numstart.dat)
c_parallel_test(posix 2 ${CDIR} ${PHASTA_BINARY_DIR}/bin/phastaC.exe)
set(compareArgs
  ${CDIR}/2-procs_case-Posix/
  ${CDIR}/2-procs_case-Posix_ref/
  0 1e-6)
c_parallel_test(compare-posix 2 ${CDIR}
  ${PHASTA_BINARY_DIR}/bin/checkphasta ${compareArgs})
c_parallel_test(restart-posix 2 ${CDIR} ${PHASTA_BINARY_DIR}/bin/phastaC.exe)
c_parallel_test(compareRestart-posix 2 ${CDIR}
  ${PHASTA_BINARY_DIR}/bin/checkphasta ${compareArgs})
c_serial_test(unlinkProcsDir-compressible ${CDIR} rm ${CDIR}/2-procs_case)
