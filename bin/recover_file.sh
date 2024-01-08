mkdir fixed
for file in $(ls *.gz); do
    #Test if the compressed file was corrupted
    if gzip -t $file; then
        echo "$file is ok"
    else 
        echo "$file is corrupt"
        echo "starting recover..."
        #get the basename
        basename=${file%.gz}
        #recover
        gzrecover $file \
        -o $basename
        gzip -f $basename > fixed/${basename}.gz
        echo "done"
    fi
done

for file in $(ls *.gz); do
    #Test if the compressed file was corrupted
    if gzip -t $file; then
        echo "$file is ok"
    else 
        echo "$file is corrupt"
    fi
done