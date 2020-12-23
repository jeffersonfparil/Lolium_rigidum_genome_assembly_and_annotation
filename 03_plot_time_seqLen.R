##########################################
### read length and a function of time ###
##########################################
args = commandArgs(trailingOnly=TRUE)
# args = c("lolium1.length.time")
fname_input = args[1] ### date/time and sequence read length tab-delimited file output of 02_fastq_time_seqLen.sh
### read input file and label the fields
dat = read.delim(fname_input, header=FALSE)
colnames(dat) = c("DATE", "SEQUENCE_LENGTH")
### parse the dates
DT = matrix(unlist((strsplit(as.character(dat$DATE), "T"))), nrow=nrow(dat), byrow=TRUE)
DT[,2] = gsub("Z", "", DT[,2])
dat$DATE = as.POSIXct(paste0(DT[,1], " ", DT[,2]))
### sort by time
dat = dat[order(dat$DATE, decreasing=FALSE), ]
### generate time field as the time elapsed
dat$TIME = dat$DATE - dat$DATE[1]
### scatter plot and polynomial regression
png(paste0(basename(fname_input), ".png"), width=1000, height=500)
par(mfrow=c(1,2))
### scatter plot
plot(x=dat$TIME/3600,
     y=(dat$SEQUENCE_LENGTH)/1000,
     xlab="Time Elapsed (hours)",
     ylab="Read Length (kb)",
     pch=20,
     col="#41b6c4",
     main="Scatterplot"
    )
grid(lty=2, col="gray")
### polynomial regression to see the relationship-ishshshshs
n_sel = 1
r2adj = summary(lm(SEQUENCE_LENGTH ~ poly(TIME, n_sel), data=dat))$adj.r.sq
for (n in 1:10){
    mod = lm(SEQUENCE_LENGTH ~ poly(TIME, n), data=dat)
    print(summary(mod)$adj.r.sq)
    if (summary(mod)$adj.r.sq >= r2adj){
        n_sel = n
        r2adj = summary(mod)$adj.r.sq
    }
}
mod = lm(SEQUENCE_LENGTH ~ poly(TIME, n_sel), data=dat)
new_time = seq(min(dat$TIME), max(dat$TIME))
y_pred = predict(mod, newdata=data.frame(TIME=new_time))
plot(x=new_time/3600, 
     y=y_pred/1000,
     xlab="Time Elapsed (hours)",
     ylab="Read Length (kb)",
     type="l",
     lty=2,
     lwd=2,
     col="#fb8072",
     main=paste0("Polynomial regression (degree=", n_sel, ")")
    )
grid(lty=2, col="gray")
dev.off()
