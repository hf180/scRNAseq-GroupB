#comparing sizes
object.size2 <- function(o1){format(object.size(o1),units="auto",standard = "SI",digits = 3)}
x <- 199
object.size2(x)

#generating big data
data1 <- matrix(rnorm(10000000),ncol=50)
write.table(data1,file = "test_manyrows.txt",row.names = FALSE,col.names = FALSE,sep = ",")

#input as datd.frame
file1 <- "test_manyrows.txt"
gc
v_times <- rbind(v_times,microbenchmark("read_df+"=data1 <-read.table(file1,header = FALSE,sep=",",colClasses= "numeric"),
                                       times = 2,unit = "s"))
v_sizes <- c(v_sizes,"inp_df"=object.size2(data1))

#view 1 
nucleotides <- c("A","T","C","G")
seq1 <- sample(nucleotides,size = 10000, replace = TRUE)
#VIEW 2
seq1_fact <- factor(seq1)
#view 3
seq1_integ <- as.integer(seq1_fact)
