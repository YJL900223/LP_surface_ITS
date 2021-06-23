##注：以下步骤均在 shell 命令行中完成，运行程序或脚本均来自 qiime
#在 otu_table.tsv 开头添加一行“# Constructed from biom file”，以便将 otu_table.tsv 转为 qiime 可识别样式
cp otu_table.txt otu_table.tsv
sed -i '1i\# Constructed from biom file' otu_table.tsv

#otu_table.tsv 转换为 otu_table.biom
biom convert -i otu_table.tsv -o otu_table.biom --table-type="OTU table" --to-json

#OTU 注释，输出结果 otu_table_tax_assignments.txt 即为注释文件
assign_taxonomy.py -i otu.fasta -r SILVA_123_SSURef_Nr99_tax_silva.fasta -t SILVA_123_SSURef_Nr99_tax_silva.tax -o ./

#添加 otu 注释信息至 biom 格式的 otu 表（otu_table.biom ）的最后一列，并将列名称命名为 taxonomy
biom add-metadata -i surface.biom --observation-metadata-fp surface_tax.tsv -o surface.biom --sc-separated taxonomy --observation-header OTUID,taxonomy 
biom add-metadata -i plant.biom --observation-metadata-fp plant_tax.tsv -o plant.biom --sc-separated taxonomy --observation-header OTUID,taxonomy
biom add-metadata -i human_skin.biom --observation-metadata-fp skin_tax.tsv -o human_skin.biom --sc-separated taxonomy --observation-header OTUID,taxonomy
biom add-metadata -i human_gut.biom --observation-metadata-fp gut_tax.tsv -o human_gut.biom --sc-separated taxonomy --observation-header OTUID,taxonomy
biom add-metadata -i indoor_air.biom --observation-metadata-fp indoor_tax.tsv -o indoor_air.biom --sc-separated taxonomy --observation-header OTUID,taxonomy

#otu_table.silva.biom 转换为 otu_table.silva.tsv
biom convert -i merged.biom -o merged.txt --to-tsv --header-key taxonomy
biom add-metadata -i otutab_norm.biom --observation-metadata-fp taxonomy2.txt \
  -o otutab_norm_tax.biom --sc-separated taxonomy \
  --observation-header OTUID,taxonomy

# txt转换为biom json格式
biom convert -i otutab_norm.txt -o otutab_norm.biom --table-type="OTU table" --to-json
# 添加物种注释
biom add-metadata -i otutab_norm.biom --observation-metadata-fp taxonomy2.txt \
  -o otutab_norm_tax.biom --sc-separated taxonomy \
  --observation-header OTUID,taxonomy
# 指定输入文件、物种注释、输出文件、注释列名、属性列名



merge_otu_tables.py -i surface.biom,plant.biom,human_skin.biom,human_gut.biom,indoor_air.biom -o merged.biom


biom summarize-table -i merged.biom -o merged_stats.txt

