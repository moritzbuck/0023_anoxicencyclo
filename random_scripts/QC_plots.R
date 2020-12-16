library(ggplot2)
library(data.table)


q_hist = fread("q_hist.csv")
md = fread("/home/moritz/kadath/data/data_submit/metadata/samples_contextual_final.csv")

setkey(md, "sample_id")
counts = q_hist[q_hist$sample %in% md$sample_id, .(nb_read_pairs = sum(count)/2) , by = sample]
counts[, country := md[sample, 'geographic location (country and/or sea)'] ]

ggplot(q_hist, aes(x=q_score, y = proportion, group=q_score))+geom_boxplot(outlier.color="white")+
  geom_jitter(alpha=0.05)+facet_grid(~direction)+theme_minimal()+theme(text = element_text(size=30))+
  xlim(0,40)+ylim(0,1)+xlab("mean Q-score of reads")


ggplot(counts, aes(x=country, y=nb_read_pairs/1000000))+geom_hline(yintercept=mean(counts$nb_read_pairs)/1000000, col="red", size=3)+
  geom_boxplot()+geom_jitter(alpha=0.7)+theme_minimal()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab('Country')+ylab('Number of read pairs (in milion)')+theme(text = element_text(size=30))
