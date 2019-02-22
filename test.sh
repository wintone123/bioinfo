for txtfile in *.txt
do 
basename=${txtfile:0:4}
    if [[ ! -f $basename/$basename.sh ]]
    then 
        echo "file not exist"
    else 
        echo "file exist"
    fi
done