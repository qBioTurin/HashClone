#!/bin/bash




OUTPUT=$(zenity  --forms --title="Hashchecker workflow"  --text="Input parameters:" --separator="," --add-entry="Window size: " --add-entry="Number of experiments: "  --add-entry="Hash table size:")
case $? in
         0)
                echo " selected.";;
         1)
                echo "No file selected."
                 exit 0;
                ;;
        -1)
                echo "An unexpected error has occurred."
                exit 0;
                ;;
                
esac

IFS=', ' read -a array <<< "$OUTPUT"

W=${array[0]}
Exp=${array[1]}
Hash=${array[2]}

if [[ ("$W" -le 10)  ||  ("$W" -ge 32) ]]; 
  then   zenity --error --text="Wrong value for Window size: $W";
  exit 0;
fi
  
if [[ ("$Exp" -le 1)  ]]; 
  then   zenity --error --text="Wrong value for Number of experiment: $Exp";
  exit 0;
fi
  
if [[ ("$Hash" -le 1)   ]]; 
  then   zenity --error --text="Wrong value for Hash table size: $Hash";
  exit 0;
fi

COUNTER=0
FILE[$Exp]=0;
echo ""> Zoutput;
echo "  =========================================">> Zoutput;
echo "|                        		     Input parameters                     		           |">> Zoutput;
echo "  =========================================">>Zoutput;
echo "">> Zoutput;
echo "Window size:  $W" >> Zoutput;
echo "">> Zoutput;
echo "Number of experiments:  $Exp" >> Zoutput;
echo "">> Zoutput;
echo "Hash table size:  $Hash" >> Zoutput;
echo "">> Zoutput;

while [  $COUNTER -lt $Exp ]; do
  echo The counter is $COUNTER
  #!/bin/sh
  FILE[$COUNTER]=`zenity --height 700 --width 400 --timeout=60 --file-selection --title="Select a File for experiment $COUNTER:"`
  case $? in
         0)
                echo "\"${FILE[$COUNTER]}\" selected."
                exit 0
                ;;
         1)
                echo "No file selected."
                exit 0
                ;;
        -1)
                echo "An unexpected error has occurred."
                exit 0
                ;;
  esac
  echo "File $COUNTER: ${FILE[$COUNTER]} ">> Zoutput;
  let COUNTER=COUNTER+1 
done
echo "">> Zoutput;
echo "==========================================">>Zoutput;

FILE=`dirname $0`/Zoutput
let HEIGHT=$Exp*100+180
zenity --text-info  --height  $HEIGHT --width 500  --timeout=60  --title="Executing Hashcheker:" --filename=$FILE \
case $? in
         0)
                echo "\"${FILE[$COUNTER]}\" selected."
                exit 0
                ;;
         1)
                echo "Execution terminated by user."
                exit 0
                ;;
        -1)
                echo "An unexpected error has occurred."
                exit 0
                ;;
  esac


