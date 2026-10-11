### Figure 1B/C

library(Seurat)
library(ggplot2)
library(patchwork)

cols_fine = c('CD4T_Tn_01_GAS5' = '#e5d25b',
              'CD4T_Tn_02_SELL' = '#599014',
              'CD4T_Tn_03_LEF1' = '#cd5d1a',
              'CD4T_Tn_04_TCF7' = '#a82d06',
              'CD4T_Tn_05_IL7R' = '#4592BF',
              'CD4T_Tcm_01_CD69' = '#d38219',
              'CD4T_Tcm_02_GPR183' = '#74a764',
              'CD4T_Tcm_03_TCF7' = '#8ca2b4',
              'CD4T_Tcm_04_S100A4' = '#cbb190',
              'CD4T_Tem' = '#e7ca8d',
              'CD4T_CTL' = '#b03e69',
              'CD4T_Treg' = '#9d9ec3',
              'CD8T_Tn_01_IL7R' = '#593202',
              'CD8T_Tn_02_TCF7' = '#467599',
              'CD8T_Tn_03_LEF1' = '#d36135',
              'CD8T_Tcm' = '#ebafa4',
              'CD8T_Tem_01_CX3CR1' = '#5e8a89',
              'CD8T_Tem_02_GZMK' = '#ecd577',
              'CD8T_Teff_01_GZMH' = '#7c606c',
              'CD8T_Teff_02_KLRG1' = '#5c6489',
              'CD8T_Tex' = '#602a2a',
              'NK_CD56_01_CCL4' = '#6fadd3',
              'NK_CD56_02_XCL1' = '#1c909a',
              'NK_CD16_01_CD44' = '#fe7e02',
              'NK_CD16_02_PRF1' = '#e4191c',
              'NK_CD16_03_IGFBP7' = '#397eb7',
              'NK_CD16_04_KLRC2' = '#4D377A',
              'NK_CD16_05_CX3CR1' = '#E5E7F3',
              'NK_CD16_06_CMC1' = '#68C3A6',
              'Bn_01_TCL1A' = '#6a7fe1',
              'Bn_02_IL4R' = '#c5aecb',
              'Bn_03_SAMHD1' = '#ca7474',
              'Bm_01_SLC25A6' = '#F1EFDE',
              'Bm_02_SAMD3' = '#8A1619',
              'Bm_03_IFITM3' = '#00589C',
              'Batyp_01_ITPR1' = '#d18e95',
              'Batyp_02_FCRLA' = '#b5b47b',
              'Bact' = '#edf2bc',
              'Bgc' = '#FFB7BE',
              'Bcyc' = '#87C48F',
              'Plasma' = '#B3DE69',
              'cMono_01_LRMDA' = '#234862',
              'cMono_02_S100A12' = '#9078d4',
              'cMono_03_TREM1' = '#A71D1D',
              'cMono_04_JUN' = '#342488',
              'ncMono_01_CALHM6' = '#797D7D',
              'ncMono_02_TBC1D8' = '#c5cc95',
              'ncMono_03_CXCR4' = '#40032b',
              'ncMono_04_BCL2' = '#dbb679',
              'ncMono_05_LYVE1' = '#c43c40',
              'cycMono' = '#4E4A8C',
              'cDC1' = '#115F5F',
              'cDC2' = '#238b5c',
              'cDC2_stress' = '#59b5b2',
              'pDC' = '#9fcccf',
              'tDC' = '#008BBD',
              'MicDC' = '#FF9F80',
              'Cycling' = '#FCCDE5')

cols_celltype = c('CD4T' = '#8DD3C7',
                  'CD8T' = '#FFFFB3',
                  'NK' = '#FB8072',
                  'B' = '#BEBADA',
                  'Plasma' = '#B3DE69',
                  'Monocyte' = '#80B1D3',
                  'DC' = '#FDB462',
                  'Cycling' = '#FCCDE5')

seurat_obj = readRDS("AD_all_annot_final.rds")

pdf(file = "Figure 1BC.pdf", width = 13, height = 6)
DimPlot(seurat_obj, group.by="fine", cols = cols_fine ,raster = T, split.by = "origin")+NoAxes()
DimPlot(seurat_obj, group.by="celltype", cols = cols_celltype, raster = T, split.by = "origin")+NoAxes()
dev.off()

### Figure S1A

g = rev(c("CD3D", "CD3E","CD4", "IL7R", "SELL", "LEF1", "GPR183", "S100A4", "GZMA", "GZMK", "KLRC2", "PRF1", "EOMES", "FOXP3", "CD8A", "CCL5", "CX3CR1", "GZMH", "FGFBP2", "GZMB", "TOX", "GNLY", "NKG7", "NCAM1", "XCL1", "CD19", "MS4A1", "TCL1A", "CRIP1", "CRIP2", "FCRL5", "JCHAIN", "CD14", "FCGR3A", "APOE", "TREM2", "CST3", "FCER1A", "CD1C", "CLEC9A", "MKI67"))

DotPlot(object = seurat_obj, features = g, col.min = -2, col.max = 2) + theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.25)) +xlab("") + ylab("") + coord_flip()+ scale_color_gradientn(colours = c("lightblue3", "lightblue", "white", "red", "red4")) + theme(axis.line=element_blank())
ggsave("Figure_S1A.pdf", width = 12, height = 9)

