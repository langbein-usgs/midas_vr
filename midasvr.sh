#!/bin/sh

if [ "$#" -lt 1 ]
then
   echo "Usage:   midasvr.sh   -f tseries.tenv | -f freeform.tenu  [-s step_file] [ -A -M -O -W value] "
   echo "          where input file tseries.tenv is in tenv2 format"
   echo "       or: freeform.tenu is in free format: station, t, e, n, u"
   echo "Purpose: Compute modified Theil-Sen estimate of velocity."
   echo "         - sample slopes between points seperated by precisely 1 year"
   echo "         Wrapper around midas_vr"
   echo "Options:  -f  name of data file ending it either .tenu or .tenv"
   echo "          -s  use step_file"
   echo "The following control computation of AVR (or VR)"
   echo "          -A  Compute AVR or VR; VR is default"
   echo "          -M  Evaluate rate variance with MAD (default) or sum of squares"
   echo "          -W  Weighting used to fit function to (A)VR data; default is 0.5"
   echo "          -O  output avr computed value in 9 files, 3 per component"
   echo "               avr_comp.dat -- computed value of A)VR  from data"
   echo "               avr_comp.pred -- predicted valus of A)VR  "
   echo "               avr_comp.fit -- fit of various models of (A)VR used to fit "
   echo "                               to (A)VR determined from the data"
   echo "Step file (free) format example:"
   echo "  RENO 2000.0000 (may have extra columns which will be copied)"
   echo "  RENO 2004.0000"
   echo "  ZUMA 2004.0000"
   echo "Output Files:"
   echo "tseries.vel   single line record with velocity solution"
   echo "tseries.renv  Residuals in tenv2 format, if input is .tenv"
   echo "tseries.renu  Residuals in free format, if input is .tenu"
   echo "tseries.step  steps matching station STID within span of data"
   exit 
fi

# set up defaults
AVRy=n
MADy=y
Scale=0.5
OUTy=n

while getopts f:s:W:AMO option
do
  case "$option"
  in
     f) infile=$OPTARG;;
     s) stepfile=$OPTARG;;
     A) AVRy=y;;
     M) MADy=n;;
     O) OUTy=y;;
     W) Scale=$OPTARG;;
     \?)  echo "Incomplete set of arguements; type midas02.sh without arguments to get documentation"
     exit 1;;
     esac
done

rm -rf MIDAS.STEPIN MIDAS.STEPOUT MIDAS.ERR MIDAS.TENV MIDAS.RENV MIDAS.TENU MIDAS.RENU MIDAS.VEL

suf=''
suf=`echo $infile | sed 's/\./ /g' | awk '{print $3}'`
if [ -z $suf ]
then 
suf=`echo $infile | sed 's/\./ /g' | awk '{print $2}'`
fi

if [ "$suf" = tenu ]
then
   ln -s $infile MIDAS.TENU
elif [ "$suf" = TENU ]
then 
   ln -s $infile MIDAS.TENU
elif [ "$suf" = TENV ]
then 
  ln -s $infile MIDAS.TENV
elif [ "$suf" = tenv  ]
then 
  ln -s $infile MIDAS.TENV
else
  echo Invalid file type, $suf
  exit
fi
sta=`echo $infile | sed 's/\./ /' | awk '{print $1}'`
if [ -f "$stepfile" ]
then
  grep $sta $stepfile > MIDAS.STEPIN
fi

rm -f avr.config
cat >> avr.config <<EOF
$MADy MADy
$AVRy AVRy
$Scale Scale for weighting
$OUTy  Output various files
EOF

# Run the code!
~/proglib/MIDAS_VR/midas_vr > MIDAS.VEL

mv MIDAS.VEL "$sta".vel
cat "$sta".vel
if [ -f  MIDAS.RENU ]
then
  mv MIDAS.RENU $sta.rneu
fi
if [ -f  MIDAS.RENV ]
then
 mv MIDAS.RENV $sta.rnev
fi
if [ -f MIDAS.STEPOUT ]
then
  mv MIDAS.STEPOUT "$sta".stepout
fi
rm -f avr.config
exit
