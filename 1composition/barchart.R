if (!require(reshape2)) install.packages("reshape2")
library(ggplot2)
library(reshape2)
data <- read.table("/Users/coellearth/Desktop/Mammary Gland Diet Project/1composition/slides.txt", header = TRUE, sep = "\t")
data_melt <- melt(data, id.vars = "CellType", variable.name = "Condition", value.name = "Value")
my_colors <- c("cornflowerblue", "firebrick")
plot <- ggplot(data_melt, aes(x = CellType, y = Value, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Mammary Gland Cell Type Composition Ration across Conditions", x = "Cell Type", y = "Composition Ratio") +
  scale_fill_manual(values = my_colors) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5))

plot + geom_text(aes(label = round(Value, 3)), 
              position = position_dodge(width = 0.9), 
              vjust = -0.5, 
              size = 2.2,
              fontface="bold")
