# qpcr
qPCR statistics

Input: 
- CSV file. Endogenous must be the first row, and columns are the treatments in 3 replicates. Change the input file name on line 15 (Template: anova_r.csv).

Statistics:
- Delta Ct (DCT) calculation by the formula DCTvalue = 2^(endogenous_control - gene_test);
- ANOVA;
- Tukey HSD.

Outputs:
- calculoqpcr.csv;
- barplot figure.


![barplot](http://amos.esalq.usp.br/hugo/barplot1.png)
Figure 1. 

\* in title: p-value < 0.05 (ANOVA)

Equal letters: differences among treatments were statistically not significant by Tukey HSD. 
