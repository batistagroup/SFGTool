####################################
### clean calculated results
####################################

for dir in ppp/* ssp sps psp ; do
    rm ${dir}/*line ${dir}/*sfg ${dir}/*dat
done
