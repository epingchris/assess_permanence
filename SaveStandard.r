SaveStandard = function(outlist, out_path, file_prefix) {
    common_var = outlist$common_var

    #summarise time series results
    summ_additionality = SummariseSim(outlist$sim_additionality)
    summ_credit = SummariseSim(outlist$sim_credit)
    summ_release = SummariseSim(outlist$sim_release[1:common_var$H, ])
    summ_aomega = SummariseSim(outlist$sim_aomega)
    summ_ep = SummariseSim(outlist$sim_ep)

    #per-year reversal risk (proportion of repetitions without positive credits)
    summ_risk = data.frame(year = 1:common_var$H, risk = apply(outlist$sim_reversal, 1, sum) / common_var$n_rep)
    summ_risk$risk[1:common_var$warmup] = NA

    #gather output objects and save
    summ_simulation = list(additionality = summ_additionality,
                           credit = summ_credit,
                           release = summ_release,
                           aomega = summ_aomega,
                           ep = summ_ep,
                           risk = summ_risk)
    summ_complete = c(common_var, summ_simulation)

    saveRDS(summ_complete, file = paste0(out_path, file_prefix, "_output.rds"))

    return(summ_complete)
}