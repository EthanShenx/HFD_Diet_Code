setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/04diffCommEnrich")

library(CellChat)

data("CellChatDB.mouse")
db.mouse <- CellChatDB.mouse

write.table(db.mouse$interaction,
          file = "CellChatDB_mouse_interaction.txt",
          row.names = FALSE, quote = FALSE, sep = "\t")
