# OXT-OXTR

### 文献 1：OXTR overexpression leads to abnormal mammary gland development in mice (*Journal of Endocrinology*, 2018)

- ++Oxtr 雌鼠（全身过表达）在未孕/早孕期即出现 **导管膨大＋过早分泌分化**，并在哺乳期迅速早衰──仔鼠因缺乳死亡 
- 机制：OXTR↑ → **催乳素-STAT5 磷酸化增强**、黄体酮↓ → 推前导管-腺泡发育时序
- 外源黄体酮可逆转提前发育，催乳素补充则可救援泌乳缺陷，说明 “OXTR-PRL/p-STAT5” 轴是核心开关

**Method**

- β-actin-OXTR 转基因鼠
- HE/免疫组化-pSTAT5
- 乳腺移植 + 激素干预 (PRL、P4) 功能验证

------

### 文献 2：Oxytocin receptor induces mammary tumorigenesis through prolactin/p-STAT5 pathway (*Cell Death & Disease*, 2021)

- 与文献 1 同一 ++Oxtr 模型，长期随访显示 **导管-腺体过度增生 → ERBB2^+ 乳腺癌**，并伴异常乳汁生成
- OXTR 透过 **PRL/p-STAT5 高激活** 驱动增殖；溴隐亭阻断 PRL 能抑瘤，提示药物可靶向这一轴心
- 揭示 OXTR 既是导管发育调控因子，也可能是 HER2 型乳癌促癌基因 

**Method**

- ++Oxtr vs WT 全程观察；乳腺全景染色评估导管/TEB 扩张
- RNA-seq & IHC 追踪 PRL/p-STAT5/ERBB2 信号

------

### 文献 3：Expression and immunolocalization of the oxytocin receptor in human lactating and non-lactating mammary glands (*Human Reproduction*, 1998)

- 在人类及狨猴乳腺中，**OXTR 不仅限于肌上皮**，在 **导管/腺体上皮** 亦检测到低水平表达，无论泌乳与否 
- 提示乳腺导管上皮本身可能是 OXT 靶点，暗示其在非泌乳条件下亦具发育或稳态功能 

**Method**

- RT-PCR 定量 OXTR mRNA

------

我们的数据：

DEG结果显示，在Stroma,Adipo,Immune cell中Oxtr显著下调

但是细胞通讯分析(cellphonedb, Nichenet, liana)都没有富集到Oxt-Oxtr这条通路

原因？

### OXT 本身极低

- 在大多数外周组织 scRNA-seq 里，几乎没有配体细胞满足 10 % 规则。
- **解决**：
  1. **降低阈值** `min.cells = 5`, `thresh = 0.05`；
  2. **合并相关亚群** 把所有少量表达 Oxt 的细胞合并成 “neuroendocrine” 伪群；
  3. **跨组织整合** 若数据只来自乳腺显然找不到 Oxt, 可把下丘脑或垂体 scRNA-seq 合并后再跑 LR。

### OXTR 低表达

- 若受体也低于阈值，则通路被整条过滤；可放宽 `min.exp = 0.05`。

