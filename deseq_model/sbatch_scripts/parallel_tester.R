#!/usr/bin/env Rscript

suppressMessages(library(BiocParallel))

print('no parallel')
register(SerialParam())
print(system.time({ bplapply(1:10, function(i) Sys.sleep(5)) }))

print('with 4')
register(MulticoreParam(workers=4))
print(system.time({ bplapply(1:10, function(i) Sys.sleep(5)) }))

print('with 10')
register(MulticoreParam(workers=10))
print(system.time({ bplapply(1:10, function(i) Sys.sleep(5)) }))
