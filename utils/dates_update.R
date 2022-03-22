
# list of dates / file names for Time Lapse - example 2_25 dataset
dates_list<-dir("cumulative_daily_AF_v4_2_25/")
write.table(dates_list,quote=F,row.names=F,col.names=F,file="dates.tsv") # manually edit order. TODO: proper date ordering using date class

# batches (batches_CMT.txt): edited manually

