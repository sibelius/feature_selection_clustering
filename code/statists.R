library(data.table)


# Perform the Wilcoxon/Mann-Whitney (W/M-W) test to compare the means of external validity indices

load("results/all6.Rdata")
mu <- results[,list(ari=mean(ari),cer=mean(cer),jaccard=mean(jaccard)),by="algorithms"]

mu_L <- results[sepVal=="L",list(ari=mean(ari),cer=mean(cer),jaccard=mean(jaccard)),by="algorithms"]
mu_M <- results[sepVal=="M",list(ari=mean(ari),cer=mean(cer),jaccard=mean(jaccard)),by="algorithms"]
mu_H <- results[sepVal=="H",list(ari=mean(ari),cer=mean(cer),jaccard=mean(jaccard)),by="algorithms"]

mu <- mu[order(-ari)]

algo <- mu$algorithms

# Apply the Wilcox/Mann-Whitney test
# if pvalue < 0.05, difference of mean of samples are significant
w_ari <- c()
w_jaccard <- c()
w_cer <- c()
for(a1 in algo)
    for(a2 in algo) {
        wt <- wilcox.test(results[algorithms==a1,ari],results[algorithms==a2,ari], "g")
        w_ari <- rbind(w_ari, data.table(a1,a2,pvalue=wt$p.value))

        wt <- wilcox.test(results[algorithms==a1,jaccard],results[algorithms==a2,jaccard], "g")
        w_jaccard <- rbind(w_jaccard, data.table(a1,a2,pvalue=wt$p.value))

        wt <- wilcox.test(results[algorithms==a1,cer],results[algorithms==a2,cer], "g")
        w_cer <- rbind(w_cer, data.table(a1,a2,pvalue=wt$p.value))
    }

mann <- function(data) {
    w_ari <- c()
    for(a1 in algo) {
        for(a2 in algo) {
            wt <- wilcox.test(data[algorithms==a1,ari],data[algorithms==a2,ari], "g")
            w_ari <- rbind(w_ari, data.table(a1,a2,pvalue=wt$p.value))
        }
    }    
    return(w_ari)
}

# generate a table of differences of mean
difference <- sapply(mu$ari, function(x) mu$ari - x)
difference_jaccard <- sapply(mu$jaccard, function(x) mu$jaccard - x)
difference <- sapply(mu$cer, function(x) mu$cer - x)

load("results/results_notKnown.Rdata")
hits <- nk[,numClust == clusterSelected,by="algorithms"]
hits_n <- hits[,sum(V1),by="algorithms"]
hits_perc <- hits[,sum(V1)/324,by="algorithms"]
