## PI16_TNFRSF21 (Stroma-LumProg)

### 文献1：The fibroblast-derived protein PI16 controls neuropathic pain (PNAS, 2020)

- PI16 is not made by neurons, glia, or immune cells but is mainly produced by fibroblasts surrounding the peripheral and central nervous system.

- PI16 promotes pain by increasing the permeability of the blood nerve barrier leading to increased immune cell infiltration
- There is evidence that PI16 regulates processing of the chemokine chemerin ([9](https://www.pnas.org/doi/full/10.1073/pnas.1913444117#core-collateral-r9)), cutaneous cathepsin K ([10](https://www.pnas.org/doi/full/10.1073/pnas.1913444117#core-collateral-r10)), and the matrix metalloprotease MMP2,PI16 plays a key role in chronic pain

Method:

- RNA-seq between G protein coupled receptor kinase 2 deficient mice (consistuitive pain) and WT mice --> peptidase inhibitor 16 (PI16) as a potential regulator of persistent pain (Samples were taken from DRG, female)

- Pi16 KO --> are protected from chronic pain
- immunoflourence
- Western blot

------

### 文献2：PI16 is a shear stress and inflammation-regulated inhibitor of MMP-2 (*Scientific Reports*, 2016)

- 血流 **高层流剪切应力** 使人冠状动脉内皮细胞 PI16 mRNA ↑119 倍、蛋白 ↑7 倍；TNF-α/IL-1β 可迅速下调该表达
- PI16 直接结合并 **抑制 MMP-2 活性**，从而降低内皮迁移并限制血管重塑
- 炎症条件（低 PI16）→ MMP-2 失控，提示 PI16 是“剪切-炎症”开关

**Method**

- micro-array & RT-qPCR
- PI16 腺病毒过表达 / siRNA 敲降 → 明胶酶谱测 MMP-2
- 人冠脉标本免疫组化定位 PI16

------

### 文献 3**Cross-tissue human fibroblast atlas reveals myofibroblast subtypes with distinct roles in immune modulation**
 *Cancer Cell*, 2024 

### **发现要点**

- - 517 份人类样本 × 269 899 个成纤维细胞单细胞 RNA-seq，定出了 20 个谱系亚型；其中 **PI16⁺“静息-储备”亚型**在 11 种组织均可见。
  - CellChat/NicheNet 预测显示：PI16⁺ 成纤维细胞与 **CX3CR1⁺ Temra/Tpex 细胞**、M2-样巨噬细胞之间的高频配体-受体通路中，**TNFRSF21（DR6） 被列为前 10% 受体靶点**，提示潜在 PI16→DR6 旁分泌信号。
  - 空间转录-多重免疫荧光证实 PI16⁺ 细胞带状分布于肿瘤边缘，而 DR6⁺ 免疫簇则在相邻免疫浸润区富集，形成“并排”微生态。
- **Method**
  - 10× Genomics 单细胞 + Harmony 跨样本整合；
  - Propeller/PROGENy 计算通路活性，CellChat & NicheNet 做配体-受体推断；
  - MIBI-TOF & CODEX 空间蛋白组验证局部共存。

------



缺少湿实验数据证明Pi16与TNFRSF21存在调控关系

我们的数据：该通讯通路在HFD条件下下调明显，但是DGE的结果显示：

- 尽管Stroma表达的Pi16有所减少(-0.3)但是未达到阈值，此外除Stroma以外，还有Adipo,Basal,Lumprog,Endo等表达Pi16，但都为不显著上调
- TNFRSF21在lumprog的表达也是非显著下调，且在7中细胞类型中都有表达



