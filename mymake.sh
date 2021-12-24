#!/bin/bash
#=======================================================================
# mymake.sh example --Zheng Gong, 2020-03-04(yy/mm/dd)
#=======================================================================

# The below line is needed to be modified if necessary.
SRC="./src"
CompilingLog="CompilationLog.txt"

#-----------------------------------------------------------------------
# Normally no need to change anything below.
PathCurrent="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
CompilingLog=$PathCurrent/$CompilingLog
TimeString=$(date  "+%Y-%m-%d %H:%M:%S")
rm -rf $CompilingLog; touch $CompilingLog
echo                                                                  | tee -a $CompilingLog
echo "!=======================*- ParaTC -*=======================!"   | tee -a $CompilingLog
echo "!                                                          !"   | tee -a $CompilingLog
echo "!     ParaTC:  Parallel Turbidity Current flow solver      !"   | tee -a $CompilingLog
echo "!     Version: 1.0                                         !"   | tee -a $CompilingLog
echo "!     Author:  Zheng Gong                                  !"   | tee -a $CompilingLog
echo "!     E-mail:  gongzheng_justin@outlook.com                !"   | tee -a $CompilingLog
echo "!                                                          !"   | tee -a $CompilingLog
echo "!====================*- Fortran 95/03 -*===================!"   | tee -a $CompilingLog
echo                                                                  | tee -a $CompilingLog
echo "  Source  Path:   "$SRC                                         | tee -a $CompilingLog
echo "  Current Path:   "$PathCurrent                                 | tee -a $CompilingLog
echo "  Compiling Time: "$TimeString                                  | tee -a $CompilingLog
echo                                                                  | tee -a $CompilingLog

EXE="ParaTC"
echo "  "$EXE"  will be compiled"                                     | tee -a $CompilingLog
echo                                                                  | tee -a $CompilingLog

echo "  Which compiler do you use? "                                  | tee -a $CompilingLog
echo "     1: Intel MPI (mpiifort)"                                   | tee -a $CompilingLog
echo "     2: gcc MPI   (mpif90). Default"                            | tee -a $CompilingLog
read -p "  Please type a compiler index (1 or 2): " id_cmp
echo    "  Please type a compiler index (1 or 2): "$id_cmp >> $CompilingLog
if [ "$id_cmp" == 1 ]; then
  CMP="intel_MPI"
else
  CMP="gcc_MPI"
fi
echo "  "$CMP"  will be used"                                         | tee -a $CompilingLog
echo                                                                  | tee -a $CompilingLog

# Compile ThirdParty
echo "  Do you want to recompile ThirdParty? "                        | tee -a $CompilingLog
echo "     1: No need, use the old ThirdParty compilation. Default"   | tee -a $CompilingLog
echo "     2: Yes, recompile ThirdParty (Recommended for first use)"  | tee -a $CompilingLog
read -p "  Please type a choice (1 or 2): " Third_flag
echo    "  Please type a choice (1 or 2): "$Third_flag >> $CompilingLog
if [ "$Third_flag" == 2 ]; then
  echo                                                                | tee -a $CompilingLog
  cd $SRC/ThirdParty
  echo "Compiling ThirdParty begins."                                 | tee -a $CompilingLog
  chmod a+x ./install_thirdParty.sh
  ./install_thirdParty.sh
  echo "Compiling ThirdParty done !"                                  | tee -a $CompilingLog
  echo                                                                | tee -a $CompilingLog
  cd ../..
else
  echo "  Choose to use the old ThirdParty compilation"               | tee -a $CompilingLog
fi
echo                                                                  | tee -a $CompilingLog

echo "  Do you want to delete temporary compiling files? "            | tee -a $CompilingLog
echo "     1: No, save them. "                                        | tee -a $CompilingLog
echo "     2: Yes,delete them. Default"                               | tee -a $CompilingLog
read -p "  Please type a choice (1 or 2): " DeleteFlag
echo    "  Please type a choice (1 or 2): "$DeleteFlag >> $CompilingLog
if [ "$DeleteFlag" != 1 ]; then
  echo "  Choose to DELETE temporary compiling files"                 | tee -a $CompilingLog
else
  echo "  Choose to SAVE temporary compiling files"                   | tee -a $CompilingLog
fi
echo                                                                  | tee -a $CompilingLog

echo  "!==================*- Compiling begins -*=================!"   | tee -a $CompilingLog
echo                                                                  | tee -a $CompilingLog
rm -fr $EXE
cd $SRC
if [ "$DeleteFlag" != 1 ]; then
  make clean >&/dev/null 
fi
make CMP=$CMP exeName=$EXE  2>&1                       | tee -a $CompilingLog
echo                                                                  | tee -a $CompilingLog
mv $EXE $PathCurrent                              
if [ $? -ne 0 ]; then
  if [ "$DeleteFlag" != 1 ]; then
    make clean >&/dev/null
  fi
  echo  $EXE" CANNOT be compiled correctly, please check !!!"         | tee -a $CompilingLog
else
  if [ "$DeleteFlag" != 1 ]; then
    make clean >&/dev/null
  fi
  echo  $EXE" has been compiled normally. Enjoy !!!"                  | tee -a $CompilingLog
  cd ..
  chmod a+x ./$EXE
fi
echo                                                                  | tee -a $CompilingLog
echo  "!===================*- Compiling ends -*==================!"   | tee -a $CompilingLog
echo                                                                  | tee -a $CompilingLog
