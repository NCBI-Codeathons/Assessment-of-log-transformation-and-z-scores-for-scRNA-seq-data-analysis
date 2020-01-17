library(Linnorm)
library(Matrix)

args <- commandArgs(trailingOnly=TRUE)
input_file <- args[1]
output_prefix <- args[2]

# Read the data
raw_data <- readMM(input_file)

linnorm_with_dot_norm <- Linnorm.Norm(raw_data, minNonZeroPortion=1e-6)
linnorm_no_dot_norm <- Linnorm(raw_data, minNonZeroPortion=1e-6)

print("Data has been normalized, now saving as mtx")

#writeMM(as(linnorm_with_dot_norm,"dgCMatrix"), file = paste0(output_prefix,"_linnorm_with_dot_norm.mtx"))
#writeMM(as(linnorm_no_dot_norm,"dgCMatrix"), file = paste0(output_prefix,"_linnorm.mtx"))

writeMM(as(linnorm_with_dot_norm,"dgCMatrix"),paste0(output_prefix,"_linnorm_with_dot_norm.mtx"))
writeMM(as(linnorm_no_dot_norm,"dgCMatrix"),paste0(output_prefix,"_linnorm.mtx"))
