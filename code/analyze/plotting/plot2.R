##======
# plot: number of genes per introgressed region
##======

print('plotting number of genes per introgressed region')

#----
# for each tag, make histogram - number of regions vs number of genes
#----

for (tag in tags) {
    print(tag)
    dir.create(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/',sep=''),showWarnings=FALSE)


    a = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plot_number_genes_by_region_',tag,'.txt',sep=''), sep='\t')
    names(a) = c('region_id', 'strain', 'chromosome', 'num_genes')

    plot = ggplot(a, aes(num_genes, fill='a')) + 
    	      geom_histogram(binwidth=1) +
	      scale_fill_viridis(discrete=TRUE) +
	      guides(fill=FALSE) +
	      scale_x_continuous(expand = c(0,0), limits=c(-1,max(a$num_genes))) +
	      xlab('number of genes per region') + ylab('number of regions') +
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(colour="black"), 
	      axis.text.y = element_text(colour="black"))

    plot + scale_y_continuous(expand = c(0,0),
               limits=c(0,max(ggplot_build(plot)$data[[1]]$count)*1.1))
    
    ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/number_genes_in_region_hist_',tag,'.pdf',sep=''), width = 12, height = 7)
}

#----
# plot average number of genes per region for all tags (or boxplot?)
#----

a = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/plot_number_genes_by_region.txt',sep=''), sep='\t')
names(a) = c('tag', 'improvement_frac', 'threshold', 'expected_length', 'expected_frac', 'mean', 'bs_lower', 'bs_upper', 'median', 'min', 'max')

ggplot(a, aes(x=tag, y=mean, fill=as.factor(improvement_frac))) + 
    	      geom_bar(stat='identity', position='dodge') +
	      ylab('average number of genes per region') + xlab('parameter set') +
	      scale_fill_viridis(discrete=TRUE) +
	      labs(fill='improvement fraction') +
	      scale_y_continuous(expand = c(0,0), limits=c(0,max(a$mean))) +
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(angle = 45,vjust = 1,hjust=1,colour="black"), 
	      axis.text.y = element_text(colour="black"))

ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/number_genes_in_region.pdf', width = 12, height = 7)


##======
# plot: lengths of all introgressed regions
##======

print('plotting lengths of introgressed regions')

#----
# for each tag, make histogram - regions vs lengths
#----

for (tag in tags) {
    print(tag)
    dir.create(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/',sep=''),showWarnings=FALSE)


    a = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plot_region_lengths_',tag,'.txt',sep=''), sep='\t')
    names(a) = c('region_id', 'strain', 'chromosome', 'region_length')

    plot = ggplot(a, aes(region_length/1000, fill='a')) + 
    	      geom_histogram(binwidth=1) +
	      xlab('region length (kb)') + ylab('number of regions') +
	      scale_fill_viridis(discrete=TRUE) +
	      guides(fill=FALSE) +
	      scale_x_continuous(expand = c(0,0), limits=c(-1,max(a$region_length)/1000)) +
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(colour="black"), 
	      axis.text.y = element_text(colour="black"))

    plot + scale_y_continuous(expand = c(0,0),
               limits=c(0,max(ggplot_build(plot)$data[[1]]$count)*1.1))
    
    ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/region_length_hist_',tag,'.pdf',sep=''), width = 12, height = 7)

}

#----
# plot average region length for all tags (or boxplot?)
#----

a = read.table(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/plot_region_lengths.txt',sep=''), sep='\t')
names(a) = c('tag', 'improvement_frac', 'threshold', 'expected_length', 'expected_frac', 'mean', 'bs_lower', 'bs_upper', 'median', 'min', 'max')

#d = dcast(d, threshold + expected_length + expected_frac ~ improvement_frac, value.var='mean')

# group by threshold
ggplot(a, aes(x=tag, y=mean/1000)) + 
    	      geom_bar(aes(fill=as.factor(threshold)), stat='identity',position='dodge') +
	      ylab('average region length (kb)') + xlab('parameter set') +
	      scale_fill_viridis(discrete=TRUE) +
	      labs(fill='improvement fraction') +
	      scale_y_continuous(expand = c(0,0), limits=c(0,max(a$mean)/1000)) +
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(angle = 45,vjust = 1,hjust=1,colour="black"), 
	      axis.text.y = element_text(colour="black"))
ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/region_length_t.pdf', width = 12, height = 7)


##======
# plot: number of introgressed bases for each strain
##======

print('plotting number of introgressed bases per strain')

	
a = read.table('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/plot_number_introgressed_bases_by_strain.txt', sep='\t')
names(a) = c('tag', 'improvement_frac', 'threshold', 'expected_length', 'expected_frac', 'strain', 'number_bases')

d = data.frame(tag = unique(a$tag), number_bases=NA, improvement_frac=NA)

#----
# for each tag, plot number of bases for each strain
#----

for (tag in tags) {
    print(tag)
    dir.create(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/',sep=''),showWarnings=FALSE)

    b = a[which(a$tag == tag),]
    d[which(d$tag == tag),]$number_bases = mean(b$number_bases)
    d[which(d$tag == tag),]$improvement_frac = b$improvement_frac[[1]]

    ggplot(b, aes(x=reorder(strain, -number_bases), y=number_bases, fill='x')) + 
    	      geom_bar(stat='identity',position='dodge') + 
	      xlab('strain') + ylab('number of bases introgressed') +
	      scale_fill_viridis(discrete=TRUE) +
	      guides(fill=FALSE) +
	      scale_y_continuous(expand = c(0,0), limits=c(0,max(b$number_bases))) +
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(angle = 45,vjust = 1,hjust=1,colour="black"), 
	      axis.text.y = element_text(colour="black"))
    
    ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/number_bases_by_strain_',tag,'.pdf',sep=''), width = 12, height = 7)

}


#----
# plot average number of bases for all tags
#----

ggplot(d, aes(x=tag, y=number_bases, fill=as.factor(improvement_frac))) + 
    	      geom_bar(stat='identity',position='dodge') + 
	      xlab('tag') + ylab('number of bases introgressed') +
	      scale_y_continuous(expand = c(0,0), limits=c(0,max(d$number_bases))) +
	      scale_fill_viridis(discrete=TRUE) +
	      labs(fill='improvement fraction') +
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(angle = 45,vjust = 1,hjust=1,colour="black"), 
	      axis.text.y = element_text(colour="black"))
    
ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/average_bases_per_strain.pdf', width = 12, height = 7)


##======
# plot: number of introgressed genes for each strain
##======

print('plotting number of introgressed genes per strain')
	
a = read.table('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/plot_number_introgressed_genes_by_strain.txt', sep='\t')
names(a) = c('tag', 'improvement_frac', 'threshold', 'expected_length', 'expected_frac', 'strain', 'number_genes')

d = data.frame(tag = unique(a$tag), number_genes=NA, improvement_frac=NA)

#----
# for each tag, plot number of genes for each strain
#----

for (tag in tags) {
    print(tag)
    dir.create(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/',sep=''),showWarnings=FALSE)

    b = a[which(a$tag == tag),]
    d[which(d$tag == tag),]$number_genes = mean(b$number_genes)
    d[which(d$tag == tag),]$improvement_frac = b$improvement_frac[[1]]
    ggplot(b, aes(x=reorder(strain, -number_genes), y=number_genes, fill='a')) + 
    	      geom_bar(stat='identity',position='dodge') + 
	      xlab('strain') + ylab('number of genes introgressed') +
	      scale_y_continuous(expand = c(0,0), limits=c(0,max(b$number_genes))) +
	      scale_fill_viridis(discrete=TRUE) +
	      guides(fill=FALSE) +
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(angle = 45,vjust = 1,hjust=1,colour="black"), 
	      axis.text.y = element_text(colour="black"))
    
    ggsave(paste('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/',tag,'/plots/number_genes_by_strain_',tag,'.pdf',sep=''), width = 12, height = 7)

}


#----
# plot average number of genes for all tags
#----

ggplot(d, aes(x=tag, y=number_genes, fill=as.factor(d$improvement_frac))) + 
    	      geom_bar(stat='identity',position='dodge') + 
	      xlab('tag') + ylab('average number of bases introgressed per strain') +
	      scale_y_continuous(expand = c(0,0), limits=c(0,max(d$number_genes))) +
	      scale_fill_viridis(discrete=TRUE) +
	      labs(fill='improvement fraction') +
              theme(panel.background=element_rect(fill="white"),
              panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
              axis.line=element_line(),
              axis.text.x = element_text(angle = 45,vjust = 1,hjust=1,colour="black"), 
	      axis.text.y = element_text(colour="black"))
    
ggsave('/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/average_genes_per_strain.pdf', width = 12, height = 7)
