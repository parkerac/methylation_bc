con <- file("log/analysis.log")
sink(con, append=TRUE)
# sink(con, append=TRUE, type="message")

# This will echo all input and not truncate 150+ character lines...
source("analysis.r", echo=TRUE, max.deparse.length=10000)

# Restore output to console
sink() 
# sink(type="message")

# And look at the log...
cat(readLines("log/analysis.log"), sep="\n")
