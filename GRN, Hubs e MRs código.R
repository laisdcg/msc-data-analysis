# ================================
# 1. Pacotes necessários
# ================================
packages <- c("tidyverse", "igraph", "ggraph", "RColorBrewer")
lapply(packages, function(pkg) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg, dependencies = TRUE)
  library(pkg, character.only = TRUE)
})
install.packages("igraph")
packageVersion("igraph")
library(igraph)
library(tidyverse)
library(ggraph)
library(RColorBrewer)
# ================================
# 2. Importar os dados
# ================================
df <- read.csv("Table of Transcription Factors.csv", stringsAsFactors = FALSE)

# Renomear colunas para facilitar
colnames(df) <- c("Gene", "Transcription Factors", "Occurrences", "Possible biological implications", "Datasets")

# Separar nome do TF (remover o motivo entre parênteses)
df$`Transcription Factors` <- gsub(" \\(.*\\)", "", df$`Transcription Factors`)

# ================================
# 3. Construir rede bipartida (Gene-TF)
# ================================
edges <- df %>% select(from = `Transcription Factors`, to = Gene, weight = Occurrences)

# Grafo direcionado TF -> Gene
g <- graph_from_data_frame(edges, directed = TRUE)

# ================================
# 4. Análise de centralidade (para encontrar Hubs)
# ================================
# Grau de entrada e saída
V(g)$in_degree <- igraph::degree(g, mode = "in")
V(g)$out_degree <- igraph::degree(g, mode = "out")
V(g)$betweenness <- betweenness(g)

# Identificar hubs (top 10 por grau total)
top_hubs <- sort(igraph::degree(g), decreasing = TRUE)[1:10]
print("Top 10 Hubs:")
print(top_hubs)

# ================================
# 5. Identificar Reguladores Mestres
# ================================
# Considera-se regulador mestre aquele TF com alto grau de saída
tf_nodes <- unique(df$TF)
tf_degrees <- igraph::degree(g, mode = "out")[tf_nodes]
master_regulators <- sort(tf_degrees, decreasing = TRUE)[1:10]
print("Reguladores Mestres:")
print(master_regulators)

# ================================
# 6. Visualização da Rede
# ================================
set.seed(123)
V(g)$type <- ifelse(V(g)$name %in% df$Gene, "Gene", "TF")
V(g)$color <- ifelse(V(g)$type == "TF", "#377eb8", "#e41a1c")
V(g)$size <- scales::rescale(igraph::degree(g), to = c(3, 10))

# Plot com ggraph
ggraph(g, layout = "fr") + 
  geom_edge_link(aes(width = weight), alpha = 0.3, color = "gray50") +
  geom_node_point(aes(color = type, size = size)) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_manual(values = c("TF" = "#377eb8", "Gene" = "#e41a1c")) +
  theme_void() +
  ggtitle("Rede de Regulação Gênica - TFs e Genes") +
  theme(legend.position = "bottom")

# ================================
# 7. Exportar Tabela de Resultados
# ================================
centrality_df <- data.frame(
  Nodo = V(g)$name,
  Tipo = V(g)$type,
  Grau_Total = igraph::degree(g),
  Grau_Entrada = V(g)$in_degree,
  Grau_Saida = V(g)$out_degree,
  Centralidade_Betweenness = V(g)$betweenness
)

write.csv(centrality_df, "centralidade_genica.csv", row.names = FALSE)

# ================================
# 8. Gráfico dos Top 10 Hubs
# ================================
top_hubs_df <- centrality_df %>%
  arrange(desc(Grau_Total)) %>%
  slice(1:10)

ggplot(top_hubs_df, aes(x = reorder(Nodo, Grau_Total), y = Grau_Total, fill = Tipo)) +
  geom_col() +
  coord_flip() +
  labs(title = "The 10 main Hubs in the Network", x = "Nodo", y = "Total Degree") +
  scale_fill_manual(values = c("TF" = "#377eb8", "Gene" = "#e41a1c")) +
  theme_minimal(base_size = 13)

# ================================
# 9. Gráfico dos Reguladores Mestres (Top 10 TFs por grau de saída)
# ================================
master_regulators_df <- centrality_df %>%
  filter(Tipo == "TF") %>%
  arrange(desc(Grau_Saida)) %>%
  slice(1:10)

ggplot(master_regulators_df, aes(x = reorder(Nodo, Grau_Saida), y = Grau_Saida)) +
  geom_col(fill = "#377eb8") +
  coord_flip() +
  labs(title = "The 10 main Master Regulators", x = "Transcription Factor", y = "Degree of Output") +
  theme_minimal(base_size = 13)



# Instale se necessário
# install.packages("visNetwork")


# Supondo que sejam 'Regulador' e 'GeneAlvo'
nodes <- data.frame(id = unique(c(df$TF, df$Gene)))
nodes$label <- nodes$id
nodes$group <- ifelse(nodes$id %in% df$TF, "TF", "Gene")

edges <- data.frame(from = df$TF, to = df$Gene)

# Visualização interativa
library(visNetwork)
visNetwork(nodes, edges) %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visGroups(groupname = "TF", color = "skyblue") %>%
  visGroups(groupname = "Gene", color = "salmon") %>%
  visLayout(randomSeed = 123) %>%
  visLegend()

