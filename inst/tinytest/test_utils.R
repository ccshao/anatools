expect_identical(rsplit("A_B_C_D_E", "_", 3), c("A_B_C", "D", "E"))
expect_identical(rs2("A_B_C_D_E", "_", 3), c("A_B", "C_D_E"))
