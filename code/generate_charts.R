pdf("subset.pdf")
ggplot(sub3, aes(x=reorder(algorithms, subset), subset)) + 
    geom_boxplot(aes(fill=algorithms)) + 
    xlab("Algoritmos") +
    ylab("Nº de Atritubos selecionados") +
    ggtitle("Nº de Atributos selecionados por cada Algoritmo")
dev.off()
