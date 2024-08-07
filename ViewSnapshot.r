#view evolution of a particular iteration
ViewSnapshot = function(ind) {
  snapshot = data.frame(additionality = sim_additionality[1:50, ind],
                        aomega = sim_aomega[1:50, ind],
                        release = sim_release[1:50, ind],
                        credit = sim_credit[1:50, ind],
                        PACT = sim_pact[1:50, ind],
                        rsched = apply(sim_schedule[[ind]], 1, sum),
                        buffer = sim_buffer[1:50, ind])
  View(snapshot)
}
