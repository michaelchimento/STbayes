test_that("extract_tv_ILV handles empty ILVs", {
    nbda_obj = STbayes::tutorial1_1

    for (ILV_type in c("asocILVdata", "intILVdata", "multiILVdata")){
        expect_message(extract_tv_ILV(nbda_obj, ILV_type))
        result <- extract_ILV(nbda_obj, ILV_type)
        expect_s3_class(result, "data.frame")
        expect_equal(ncol(result), 1)
        expect_equal(colnames(result), c("ILVabsent"))
        expect_equal(sum(result$ILVabsent), 0)
    }
})

test_that("extract_ILV handles asoc_ilv, int_ilv", {
    nbda_obj = STbayes::tutorial2_2

    for (ILV_type in c("asocILVdata", "intILVdata")){
        result <- extract_tv_ILV(nbda_obj, "asocILVdata")
        expect_s3_class(result, "data.frame")
        expect_equal(ncol(result), 2)
        expect_equal(colnames(result), c("treatment", "maleTV"))
    }
})

test_that("extract_ILV errors if slot name is wrong", {
    nbda_obj = STbayes::tutorial2_2
    expect_error(extract_tv_ILV(nbda_obj, "nonexistentSlot"))
})
