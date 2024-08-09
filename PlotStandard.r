#plot summary of time series for credit, EP and reversal risk
PlotStandard = function(summary_out, range_credit, out_path, show_axis_title = F, file_type = "png", file_prefix) {
    type = summary_out$type
    H = summary_out$H
    year_label = c(1, seq(H / 10, H, by = H / 10))
    t_max_x = summary_out$t_max - summary_out$t0 + 1
    
    #a. Credits
    ggplot(summary_out$credit, aes(x = year)) +
        geom_ribbon(aes(ymin = p05, ymax = p95), fill = "lightgray", alpha = 0.5) +
        geom_ribbon(aes(ymin = p25, ymax = p75), fill = "darkgray", alpha = 0.5) +
        geom_line(aes(y = median), size = 1, color = "black") +
        geom_hline(yintercept = 0, size = 0.5, color = "black") +
        {if(type == "real") geom_vline(xintercept = t_max_x, size = 0.5, color = "black", lty = "dotted")} +
        scale_x_continuous(breaks = year_label, labels = year_label) +
        scale_y_continuous(limits = range_credit) +
        labs(x = ifelse(show_axis_title, "Year", ""),
             y = ifelse(show_axis_title, expression(paste("Credits (Mg ", CO[2], ")", sep = "")), "")) +
        theme_bw() +
        theme(axis.line = element_line(linewidth = 0.5),
              panel.grid = element_blank(),
              axis.title = element_text(size = 26),
              axis.text.x = element_text(size = 18, vjust = 0.5),
              axis.text.y = element_text(size = 22),
              plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"))
    ggsave(paste0(out_path, file_prefix, "_a_credit.", file_type), width = 12, height = 9, unit = "cm")

    #b. EP
    ggplot(summary_out$ep, aes(x = year)) +
        geom_ribbon(aes(ymin = p05, ymax = p95), fill = "lightgray", alpha = 0.5) +
        geom_ribbon(aes(ymin = p25, ymax = p75), fill = "darkgray", alpha = 0.5) +
        geom_line(aes(y = median), size = 1, color = "black") +
        {if(type == "real") geom_vline(xintercept = t_max_x, size = 0.5, color = "black", lty = "dotted")} +
        scale_x_continuous(breaks = year_label, labels = year_label) +
        scale_y_continuous(limits = c(0, 1)) +
        labs(x = ifelse(show_axis_title, "Year", ""),
             y = ifelse(show_axis_title, "EP", "")) +
        theme_bw() +
        theme(axis.line = element_line(linewidth = 0.5),
              panel.grid = element_blank(),
              axis.title = element_text(size = 26),
              axis.text.x = element_text(size = 18, vjust = 0.5),
              axis.text.y = element_text(size = 22),
              plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"))
    ggsave(paste0(out_path, file_prefix, "_b_ep.", file_type), width = 12, height = 9, unit = "cm")

    #c. Reversal risk
    ggplot(summary_out$risk, aes(x = year)) +
        geom_line(aes(y = risk), size = 1) +
        geom_hline(yintercept = 0.05, color = "red", lty = "dashed", size = 0.5) +
        {if(type == "real") geom_vline(xintercept = t_max_x, size = 0.5, color = "black", lty = "dotted")} +
        scale_x_continuous(breaks = year_label, labels = year_label) +
        scale_y_continuous(limits = c(0, 0.5)) +
        labs(x = ifelse(show_axis_title, "Year", ""),
             y = ifelse(show_axis_title, "Reversal risk", "")) +
        theme_bw() +
        theme(axis.line = element_line(linewidth = 0.5),
              panel.grid = element_blank(),
              axis.title = element_text(size = 26),
              axis.text.x = element_text(size = 18, vjust = 0.5),
              axis.text.y = element_text(size = 22),
              plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"))
    ggsave(paste0(out_path, file_prefix, "_c_risk.", file_type), width = 12, height = 9, unit = "cm")
}

