# ============================================================
# NeuroLLPSAI V3.0 产品化竞赛原型
# 场景：临床科研 / 转化研究 / 分子检测结果辅助解释
# 工作流：数据导入 -> 数据质控 -> AI推理 -> 风险解释 -> 报告输出
# 仅用于科研与竞赛展示，不用于临床诊疗决策
# ============================================================

library(shiny)
library(randomForestSRC)

# ------------------------------------------------------------
# 1. 模型文件
# GitHub / Connect Cloud 根目录需包含：
# app.R
# manifest.json
# rsf_model.rds
# seven_gene_model.rds
# ------------------------------------------------------------
paper_model_file <- "seven_gene_model.rds"
rsf_model_file <- "rsf_model.rds"

if (!file.exists(paper_model_file)) stop("找不到 seven_gene_model.rds")
if (!file.exists(rsf_model_file)) stop("找不到 rsf_model.rds")

paper <- readRDS(paper_model_file)
rsf_model <- readRDS(rsf_model_file)

genes <- paper$genes
model_coef <- paper$model_coef[genes]
risk_cutoff <- as.numeric(paper$risk_cutoff)

expr_train <- paper$expr_train[genes, , drop = FALSE]
expr_valid <- paper$expr_valid[genes, , drop = FALSE]
train_risk_df <- paper$train_risk_df
valid_risk_df <- paper$valid_risk_df

gene_min <- apply(expr_train, 1, min, na.rm = TRUE)
gene_max <- apply(expr_train, 1, max, na.rm = TRUE)

rsf_train_score <- as.numeric(rsf_model$predicted)
rsf_cutoff <- median(rsf_train_score, na.rm = TRUE)

# ------------------------------------------------------------
# 2. 真实验证集中的典型演示病例
# ------------------------------------------------------------
valid_x <- as.data.frame(t(expr_valid[genes, , drop = FALSE]))
valid_x <- valid_x[, genes, drop = FALSE]
valid_samples <- rownames(valid_x)

paper_valid_match <- valid_risk_df[
  match(valid_samples, valid_risk_df$sample),
  ,
  drop = FALSE
]

rsf_valid_pred <- predict(rsf_model, newdata = valid_x)
rsf_valid_score <- as.numeric(rsf_valid_pred$predicted)

cox_valid_score <- paper_valid_match$riskScore
cox_valid_high <- cox_valid_score > risk_cutoff
rsf_valid_high <- rsf_valid_score > rsf_cutoff

within_training_range <- apply(
  valid_x,
  1,
  function(v) all(v >= gene_min[genes] & v <= gene_max[genes])
)

rsf_valid_times <- as.numeric(rsf_valid_pred$time.interest)
idx_5y <- which.min(abs(rsf_valid_times - 365 * 5))

if (is.matrix(rsf_valid_pred$survival)) {
  valid_surv_5y <- as.numeric(rsf_valid_pred$survival[, idx_5y])
} else {
  valid_surv_5y <- rep(NA_real_, nrow(valid_x))
}

low_candidates <- which(!cox_valid_high & !rsf_valid_high & within_training_range)
high_candidates <- which(cox_valid_high & rsf_valid_high & within_training_range)

pick_middle_case <- function(indices, surv5) {
  if (length(indices) == 0) return(NA_integer_)
  vals <- surv5[indices]
  target <- median(vals, na.rm = TRUE)
  indices[which.min(abs(vals - target))]
}

low_idx <- pick_middle_case(low_candidates, valid_surv_5y)
high_idx <- pick_middle_case(high_candidates, valid_surv_5y)

if (is.na(low_idx)) low_idx <- which.min(cox_valid_score)
if (is.na(high_idx)) high_idx <- which.max(cox_valid_score)

low_values <- setNames(as.numeric(expr_valid[, valid_samples[low_idx]]), genes)
high_values <- setNames(as.numeric(expr_valid[, valid_samples[high_idx]]), genes)

# ------------------------------------------------------------
# 3. 已完成并核对的模型验证指标
# ------------------------------------------------------------
internal_perf <- data.frame(
  Time = c("1年", "3年", "5年"),
  SevenGene = c(0.758, 0.871, 0.882),
  RSF = c(0.798, 0.880, 0.888),
  stringsAsFactors = FALSE
)

cindex_values <- c(
  SevenGene = 0.832,
  RSF = 0.843
)

external_auc <- c(
  "1年" = 0.843,
  "3年" = 0.773,
  "5年" = 0.787
)

if (!is.null(rsf_model$importance)) {
  vimp <- data.frame(
    Gene = names(rsf_model$importance),
    Importance = as.numeric(rsf_model$importance),
    stringsAsFactors = FALSE
  )
  vimp <- vimp[order(vimp$Importance, decreasing = TRUE), ]
} else {
  vimp <- data.frame(
    Gene = c("CHD5","NRCAM","ABLIM3","ALDOC","NTRK1","MAP1LC3A","CHEK1"),
    Importance = c(
      0.05894561, 0.04150449, 0.03901405, 0.02457812,
      0.02376363, 0.01330599, 0.01027306
    ),
    stringsAsFactors = FALSE
  )
}

qpcr_result <- c(
  "NRCAM" = "下降",
  "MAP1LC3A" = "下降",
  "NTRK1" = "下降",
  "ALDOC" = "下降",
  "CHEK1" = "升高",
  "ABLIM3" = "下降趋势",
  "CHD5" = "未可靠定量"
)

# ------------------------------------------------------------
# 4. 工具函数
# ------------------------------------------------------------
clamp <- function(x, lo = 0, hi = 100) {
  pmax(lo, pmin(hi, x))
}

percentile_of <- function(x, reference) {
  100 * mean(reference <= x, na.rm = TRUE)
}

safe_html <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

calc_prediction <- function(values) {
  values <- as.numeric(values)
  names(values) <- genes

  outside <- genes[
    values < gene_min[genes] |
    values > gene_max[genes]
  ]

  cox_score <- sum(values * model_coef[genes])
  cox_high <- cox_score > risk_cutoff
  cox_pct <- percentile_of(cox_score, train_risk_df$riskScore)

  newdata <- as.data.frame(as.list(values))
  newdata <- newdata[, genes, drop = FALSE]

  rsf_pred <- predict(rsf_model, newdata = newdata)
  rsf_score <- as.numeric(rsf_pred$predicted)[1]
  rsf_high <- rsf_score > rsf_cutoff
  rsf_pct <- percentile_of(rsf_score, rsf_train_score)

  pred_times <- as.numeric(rsf_pred$time.interest)
  surv_curve <- rsf_pred$survival
  if (!is.null(dim(surv_curve))) surv_curve <- surv_curve[1, ]
  surv_curve <- as.numeric(surv_curve)

  get_surv <- function(day) {
    idx <- which.min(abs(pred_times - day))
    max(min(surv_curve[idx], 1), 0)
  }

  list(
    values = values,
    outside = outside,
    cox_score = cox_score,
    cox_high = cox_high,
    cox_pct = cox_pct,
    rsf_score = rsf_score,
    rsf_high = rsf_high,
    rsf_pct = rsf_pct,
    pred_times = pred_times,
    surv_curve = surv_curve,
    surv_1y = get_surv(365),
    surv_3y = get_surv(365 * 3),
    surv_5y = get_surv(365 * 5)
  )
}

overall_text <- function(x) {
  agree <- identical(x$cox_high, x$rsf_high)

  if (agree && x$cox_high) {
    list(
      level = "较高风险方向",
      brief = "七基因模型与RSF均提示较高风险方向。",
      class = "high"
    )
  } else if (agree && !x$cox_high) {
    list(
      level = "较低风险方向",
      brief = "七基因模型与RSF均提示较低风险方向。",
      class = "low"
    )
  } else {
    list(
      level = "模型意见不完全一致",
      brief = "两种模型对该样本的风险方向判断不同，应结合完整临床信息谨慎解释。",
      class = "mixed"
    )
  }
}

build_qc <- function(values, source_name = "当前输入") {
  values <- as.numeric(values)
  names(values) <- genes

  finite_ok <- all(is.finite(values))
  complete_ok <- length(values) == length(genes) && all(genes %in% names(values))
  outside <- genes[
    is.finite(values) &
    (values < gene_min[genes] | values > gene_max[genes])
  ]
  range_ok <- finite_ok && length(outside) == 0

  list(
    source = source_name,
    complete_ok = complete_ok,
    finite_ok = finite_ok,
    range_ok = range_ok,
    outside = outside
  )
}

status_badge <- function(ok, good, bad, warn = FALSE) {
  if (warn) {
    div(class = "qc-badge qc-warn", "!", div(class = "qc-text", bad))
  } else if (ok) {
    div(class = "qc-badge qc-ok", "✓", div(class = "qc-text", good))
  } else {
    div(class = "qc-badge qc-bad", "×", div(class = "qc-text", bad))
  }
}

metric_card <- function(number, label, note = NULL) {
  div(
    class = "metric-card",
    div(class = "metric-number", number),
    div(class = "metric-label", label),
    if (!is.null(note)) div(class = "metric-note", note)
  )
}

risk_gauge <- function(percentile) {
  p <- clamp(percentile)
  div(
    class = "risk-gauge-wrap",
    div(
      class = "risk-gauge",
      div(class = "risk-marker", style = paste0("left:calc(", p, "% - 7px);"))
    ),
    div(
      class = "risk-label-row",
      span("较低风险位置"),
      strong(paste0("训练集第 ", round(p), " 百分位")),
      span("较高风险位置")
    )
  )
}

score_bar <- function(label, value, max_value = 1) {
  pct <- clamp(100 * value / max_value)
  div(
    class = "score-block",
    div(
      class = "score-head",
      span(label),
      strong(sprintf("%.3f", value))
    ),
    div(
      class = "bar-track",
      div(class = "bar-fill", style = paste0("width:", pct, "%;"))
    )
  )
}

# ------------------------------------------------------------
# 5. CSS
# ------------------------------------------------------------
css_text <- "
*,
*::before,
*::after{
  box-sizing:border-box;
}
html,
body{
  width:100%;
  max-width:100%;
  overflow-x:hidden;
}
body{
  background:#f4f7fb;
  color:#172033;
  font-family:'Microsoft YaHei','PingFang SC','Noto Sans CJK SC',Arial,sans-serif;
}
.row,
[class*='col-'],
.card,
.metric-card,
.flow-card,
.upload-zone,
.report-shell{
  min-width:0;
  max-width:100%;
}
.card,
.metric-card,
.flow-card,
.upload-zone,
.report-shell{
  overflow-wrap:anywhere;
}
.navbar{
  background:#102f4f !important;
  border:none !important;
  box-shadow:0 2px 14px rgba(16,47,79,.14);
}
.navbar-default .navbar-brand,
.navbar-default .navbar-nav>li>a{
  color:#fff !important;
  font-weight:650;
}
.navbar-default .navbar-brand{
  cursor:pointer;
  transition:background .15s ease,color .15s ease;
}
.navbar-default .navbar-brand:hover,
.navbar-default .navbar-brand:focus{
  background:#0c2742 !important;
  color:#fff !important;
}
.navbar-default .navbar-nav>.active>a,
.navbar-default .navbar-nav>.active>a:hover{
  background:#1d5f8d !important;
  color:#fff !important;
}
.container-fluid{
  width:100%;
  max-width:1340px;
  padding-left:20px;
  padding-right:20px;
}
.form-control,
.progress,
.shiny-input-container,
.shiny-file-input-progress{
  max-width:100%;
}
.hero{
  margin:20px 0 18px 0;
  padding:42px 40px;
  border-radius:22px;
  color:white;
  background:linear-gradient(135deg,#102f4f 0%,#1f6b98 100%);
  box-shadow:0 12px 30px rgba(16,47,79,.18);
}
.hero-kicker{
  font-size:13px;
  letter-spacing:1.6px;
  opacity:.80;
  margin-bottom:9px;
}
.hero h1{
  margin:0 0 12px 0;
  font-size:40px;
  font-weight:850;
}
.hero-sub{
  font-size:19px;
  line-height:1.65;
  margin-bottom:10px;
  max-width:900px;
}
.hero-mission{
  margin-top:4px;
  font-size:16px;
  line-height:1.7;
  font-weight:700;
  opacity:.96;
}
.hero-note{
  margin-top:5px;
  font-size:13px;
  opacity:.82;
}
.hero-actions{
  margin-top:22px;
  display:flex;
  gap:10px;
  flex-wrap:wrap;
}
.hero-btn{
  display:inline-block;
  padding:11px 18px;
  border-radius:10px;
  background:white;
  color:#174d70 !important;
  font-weight:800;
  text-decoration:none !important;
  cursor:pointer;
}
.hero-btn-secondary{
  display:inline-block;
  padding:11px 18px;
  border-radius:10px;
  border:1px solid rgba(255,255,255,.55);
  color:white !important;
  font-weight:750;
  text-decoration:none !important;
  cursor:pointer;
}
.card{
  background:#fff;
  border:1px solid #e7edf3;
  border-radius:18px;
  padding:22px;
  margin-bottom:18px;
  box-shadow:0 5px 18px rgba(28,61,91,.07);
}
.section-title{
  font-size:22px;
  font-weight:820;
  margin:0 0 14px 0;
}
.section-sub{
  color:#667585;
  line-height:1.7;
  margin-bottom:15px;
}
.metric-row{
  display:grid;
  grid-template-columns:repeat(4,1fr);
  gap:16px;
  margin-bottom:18px;
}
.metric-card{
  background:white;
  border:1px solid #e7edf3;
  border-radius:18px;
  padding:22px 16px;
  text-align:center;
  box-shadow:0 5px 18px rgba(28,61,91,.07);
}
.metric-number{
  color:#1e638f;
  font-size:34px;
  line-height:1.1;
  font-weight:850;
}
.metric-label{
  margin-top:8px;
  font-size:14px;
  color:#4f5f6f;
}
.metric-note{
  margin-top:5px;
  font-size:11px;
  color:#8a96a3;
}
.workflow{
  display:flex;
  flex-wrap:wrap;
  gap:12px;
  width:100%;
}
.flow-card{
  position:relative;
  flex:1 1 205px;
  min-width:0;
  padding:18px 14px;
  background:#f7fbfe;
  border:1px solid #dfeaf2;
  border-radius:15px;
  min-height:125px;
}
.flow-num{
  color:#1e638f;
  font-size:12px;
  font-weight:850;
  letter-spacing:1px;
}
.flow-title{
  margin:8px 0 5px 0;
  font-size:17px;
  font-weight:820;
}
.flow-text{
  font-size:12px;
  color:#667585;
  line-height:1.6;
}
.three-grid{
  display:grid;
  grid-template-columns:repeat(3,minmax(0,1fr));
  gap:16px;
}
.feature-card{
  border-radius:15px;
  padding:18px;
  background:#f8fbfd;
  border:1px solid #e1ebf2;
}
.feature-title{
  margin:4px 0 7px 0;
  font-size:18px;
  font-weight:820;
}
.feature-text{
  color:#5e6d7c;
  line-height:1.65;
  font-size:14px;
}
.audience-tag{
  display:inline-flex;
  padding:7px 11px;
  background:#edf5fb;
  border-radius:999px;
  color:#315e7e;
  font-size:13px;
  font-weight:750;
  margin:4px 5px 4px 0;
}
.mode-tabs{
  display:grid;
  grid-template-columns:repeat(2,minmax(0,1fr));
  gap:10px;
  margin-bottom:16px;
}
.mode-tabs .btn,
.mode-tabs .shiny-download-link{
  width:100%;
  min-width:0;
  min-height:40px;
  white-space:normal;
  overflow:visible;
  text-overflow:clip;
  line-height:1.25;
  cursor:pointer;
}
.mode-tabs .shiny-download-link{
  grid-column:1 / -1;
  display:flex;
  align-items:center;
  justify-content:center;
}
.mode-btn{
  padding:13px 12px;
  border:1px solid #d8e3eb;
  background:#f9fbfd;
  border-radius:12px;
  text-align:center;
  font-weight:750;
  cursor:pointer;
}
.advanced-panel{
  margin-top:14px;
  border-top:1px solid #edf1f4;
  padding-top:10px;
}
.advanced-summary{
  cursor:pointer;
  user-select:none;
  -webkit-user-select:none;
  display:flex;
  align-items:center;
  justify-content:space-between;
  padding:11px 12px;
  border-radius:10px;
  color:#27485f;
  background:#f7fafc;
  border:1px solid #e3ebf1;
  font-weight:750;
  list-style:none;
  transition:background .15s ease,border-color .15s ease;
}
.advanced-summary:hover{
  background:#eef6fc;
  border-color:#c7ddea;
}
.advanced-summary::-webkit-details-marker{
  display:none;
}
.advanced-summary::after{
  content:'▾';
  color:#1e638f;
  font-size:15px;
  margin-left:10px;
}
details[open] > .advanced-summary::after{
  content:'▴';
}
.advanced-panel[open] > .advanced-summary{
  margin-bottom:10px;
}
.upload-zone{
  padding:16px;
  background:#f7fbfe;
  border:1px dashed #aac7da;
  border-radius:14px;
}
.info-strip{
  padding:13px 15px;
  border-radius:12px;
  background:#eef6fc;
  border:1px solid #d7e8f4;
  line-height:1.65;
}
.warn-strip{
  padding:13px 15px;
  border-radius:12px;
  background:#fff7e8;
  border:1px solid #efdcae;
  line-height:1.65;
  margin-bottom:13px;
}
.good-strip{
  padding:13px 15px;
  border-radius:12px;
  background:#edf8f2;
  border:1px solid #cce7d7;
  line-height:1.65;
}
.qc-grid{
  display:grid;
  grid-template-columns:repeat(4,1fr);
  gap:10px;
  margin:13px 0 4px 0;
}
.qc-badge{
  padding:13px;
  border-radius:13px;
  font-size:22px;
  font-weight:850;
  text-align:center;
}
.qc-text{
  margin-top:5px;
  font-size:11px;
  font-weight:650;
}
.qc-ok{background:#edf8f2;color:#27734c;border:1px solid #cce7d7;}
.qc-warn{background:#fff7e8;color:#9a6b1d;border:1px solid #efdcae;}
.qc-bad{background:#fdeeee;color:#a34b4b;border:1px solid #efcccc;}
.qc-pending{
  background:#f4f6f8;
  color:#85919d;
  border:1px solid #dfe5ea;
}
.qc-pending .qc-symbol{
  display:inline-block;
  animation:qcPulse 1s ease-in-out infinite;
}
@keyframes qcPulse{
  0%,100%{opacity:.35;}
  50%{opacity:1;}
}
.data-table{
  width:100%;
  border-collapse:collapse;
  margin-top:10px;
}
.data-table th,.data-table td{
  padding:8px 10px;
  border-bottom:1px solid #edf1f4;
  text-align:left;
}
.data-table th{
  background:#f6f9fb;
  font-size:12px;
}
.result-hero{
  padding:24px;
  border-radius:17px;
  background:#f7fafc;
  border:1px solid #e1e9ef;
  margin-bottom:14px;
}
.result-eyebrow{
  font-size:13px;
  color:#637383;
  margin-bottom:4px;
}
.result-main{
  font-size:32px;
  font-weight:850;
  margin-bottom:5px;
}
.result-sub{
  color:#5d6b79;
  line-height:1.65;
}
.model-grid{
  display:grid;
  grid-template-columns:1fr 1fr;
  gap:14px;
  margin-bottom:14px;
}
.model-card{
  padding:18px;
  border-radius:15px;
  background:#fff;
  border:1px solid #e0e8ef;
}
.model-name{
  font-size:15px;
  color:#617080;
}
.model-judgement{
  font-size:25px;
  font-weight:820;
  margin:4px 0;
}
.risk-gauge-wrap{
  margin-top:11px;
}
.risk-gauge{
  position:relative;
  height:10px;
  border-radius:999px;
  background:linear-gradient(90deg,#d8edf8 0%,#eef2f5 50%,#f2d7d4 100%);
}
.risk-marker{
  position:absolute;
  top:-4px;
  width:14px;
  height:18px;
  border-radius:5px;
  background:#163f62;
  box-shadow:0 1px 5px rgba(0,0,0,.2);
}
.risk-label-row{
  display:flex;
  justify-content:space-between;
  gap:8px;
  font-size:11px;
  color:#7b8793;
  margin-top:6px;
}
.risk-label-row strong{
  color:#384757;
}
.survival-summary{
  display:grid;
  grid-template-columns:repeat(3,minmax(0,1fr));
  gap:12px;
  margin:14px 0 6px 0;
}
.surv-mini{
  padding:14px 10px;
  background:#f6f9fb;
  border-radius:13px;
  text-align:center;
}
.surv-year{
  font-size:12px;
  color:#687786;
}
.surv-value{
  font-size:24px;
  font-weight:820;
  margin-top:4px;
}
.report-shell{
  background:#fff;
  border:1px solid #dfe7ee;
  border-radius:18px;
  padding:28px;
  box-shadow:0 6px 18px rgba(28,61,91,.08);
  margin-bottom:18px;
}
.report-header{
  display:flex;
  justify-content:space-between;
  gap:20px;
  align-items:flex-start;
  border-bottom:2px solid #e7edf3;
  padding-bottom:17px;
  margin-bottom:18px;
}
.report-title{
  font-size:27px;
  font-weight:850;
}
.report-meta{
  text-align:right;
  color:#6a7886;
  font-size:12px;
  line-height:1.65;
}
.report-section{
  padding:15px 0;
  border-bottom:1px solid #edf1f4;
}
.report-section:last-child{
  border-bottom:none;
}
.report-section-title{
  font-size:17px;
  font-weight:820;
  margin-bottom:10px;
}
.report-actions{
  display:flex;
  gap:10px;
  flex-wrap:wrap;
  margin-bottom:18px;
}
.score-block{
  margin-bottom:16px;
}
.score-head{
  display:flex;
  justify-content:space-between;
  align-items:center;
  margin-bottom:6px;
}
.bar-track{
  width:100%;
  height:12px;
  background:#edf1f4;
  border-radius:999px;
  overflow:hidden;
}
.bar-fill{
  height:100%;
  border-radius:999px;
  background:linear-gradient(90deg,#7bb3d5,#1e638f);
}
.pair-row{
  display:grid;
  grid-template-columns:90px 1fr;
  gap:12px;
  align-items:start;
  padding:11px 0;
  border-bottom:1px solid #edf1f4;
}
.pair-time{
  font-weight:820;
}
.pair-values{
  display:grid;
  grid-template-columns:1fr 1fr;
  gap:10px;
}
.pair-chip{
  background:#f6f9fb;
  border-radius:11px;
  padding:10px 12px;
}
.pair-chip strong{
  float:right;
}
.vimp-row{
  margin-bottom:11px;
}
.vimp-head{
  display:flex;
  justify-content:space-between;
  margin-bottom:5px;
  font-size:13px;
}
.vimp-track{
  height:10px;
  background:#edf1f4;
  border-radius:999px;
  overflow:hidden;
}
.vimp-fill{
  height:100%;
  border-radius:999px;
  background:#1e638f;
}
.evidence-grid{
  display:grid;
  grid-template-columns:repeat(4,1fr);
  gap:14px;
}
.evidence-card{
  padding:16px;
  border-radius:14px;
  background:#f8fbfd;
  border:1px solid #e1eaf1;
}
.evidence-number{
  font-size:27px;
  font-weight:820;
  color:#1e638f;
}
.evidence-label{
  margin-top:5px;
  font-size:13px;
  color:#536372;
}
.roadmap{
  display:flex;
  flex-wrap:wrap;
  align-items:center;
  gap:9px;
}
.road-step{
  padding:12px 14px;
  border-radius:12px;
  background:#f6f9fb;
  border:1px solid #dfe7ee;
  font-weight:750;
}
.road-step-current{
  background:#e9f4fb;
  border-color:#bcd9ea;
  color:#1e638f;
}
.road-arrow{
  color:#8796a4;
}
.small-note{
  color:#74818e;
  font-size:12px;
  line-height:1.65;
}
.footer-box{
  color:#687786;
  font-size:12px;
  line-height:1.7;
}
.btn-primary{
  background:#1e638f;
  border-color:#1e638f;
}
.btn-primary:hover{
  background:#164e72;
  border-color:#164e72;
}
.form-control{
  border-radius:10px;
}


.maturity-strip{
  display:flex;
  align-items:flex-start;
  gap:14px;
  padding:16px 18px;
  margin-bottom:18px;
  border-radius:15px;
  background:linear-gradient(90deg,#edf6fb,#f8fbfd);
  border:1px solid #d4e6f1;
}
.maturity-label{
  flex:0 0 auto;
  padding:6px 10px;
  border-radius:999px;
  background:#1e638f;
  color:white;
  font-size:12px;
  font-weight:800;
}
.maturity-text{
  color:#435566;
  line-height:1.7;
  font-size:14px;
}
.tech-principle{
  display:grid;
  grid-template-columns:1fr 44px 1fr;
  gap:14px;
  align-items:stretch;
}
.tech-stage{
  background:#f8fbfd;
  border:1px solid #dfeaf2;
  border-radius:16px;
  padding:20px;
}
.tech-stage-title{
  font-size:18px;
  font-weight:820;
  margin-bottom:5px;
}
.tech-stage-sub{
  font-size:12px;
  color:#71808e;
  margin-bottom:15px;
  line-height:1.6;
}
.tech-step{
  display:grid;
  grid-template-columns:34px 1fr;
  gap:10px;
  align-items:start;
  padding:10px 0;
  border-bottom:1px solid #edf1f4;
}
.tech-step:last-child{
  border-bottom:none;
}
.tech-num{
  width:30px;
  height:30px;
  display:flex;
  align-items:center;
  justify-content:center;
  border-radius:9px;
  background:#e8f3fa;
  color:#1e638f;
  font-weight:850;
  font-size:12px;
}
.tech-step-title{
  font-weight:800;
  margin-bottom:3px;
}
.tech-step-text{
  color:#657482;
  font-size:12px;
  line-height:1.55;
}
.tech-arrow-column{
  display:flex;
  align-items:center;
  justify-content:center;
  color:#7b91a5;
  font-size:30px;
  font-weight:300;
}
.tech-output-grid{
  display:grid;
  grid-template-columns:repeat(2,minmax(0,1fr));
  gap:9px;
  margin-top:12px;
}
.tech-output{
  padding:10px 11px;
  border-radius:11px;
  background:white;
  border:1px solid #e1e9ef;
  font-size:12px;
  font-weight:750;
  color:#405363;
}
.transparency-note{
  margin-top:15px;
  padding:13px 15px;
  border-left:4px solid #1e638f;
  background:#f5f9fc;
  color:#536473;
  line-height:1.7;
  font-size:12px;
}
@media(max-width:900px){
  .tech-principle{
    grid-template-columns:1fr;
  }
  .tech-arrow-column{
    transform:rotate(90deg);
    min-height:25px;
  }
}

.site-footer{
  max-width:1340px;
  margin:8px auto 24px auto;
  padding:13px 20px;
  color:#74818e;
  font-size:12px;
  line-height:1.7;
  text-align:center;
}

@media(max-width:950px){
  .metric-row,.evidence-grid{grid-template-columns:repeat(2,1fr);}
  .three-grid,.model-grid{grid-template-columns:1fr;}
  .qc-grid{grid-template-columns:1fr 1fr;}
}
@media(max-width:600px){
  .metric-row,.evidence-grid,.survival-summary,.qc-grid{grid-template-columns:1fr;}
  .flow-card{flex-basis:100%;}
  .hero{padding:28px 24px;}
  .hero h1{font-size:30px;}
  .report-header{display:block;}
  .report-meta{text-align:left;margin-top:9px;}
}
@media print{
  .navbar,.report-actions,.shiny-notification,.no-print{display:none !important;}
  body{background:white;}
  .container-fluid{max-width:none;}
  .report-shell{box-shadow:none;border:none;}
}
"

# ------------------------------------------------------------
# 6. UI
# ------------------------------------------------------------
ui <- navbarPage(
  title = "NeuroLLPSAI",
  id = "tabs",
  header = tags$head(
    tags$style(HTML(css_text)),
    tags$script(HTML("
      $(document).on('click', '.navbar-brand', function(e){
        e.preventDefault();
        Shiny.setInputValue('brand_home', Date.now(), {priority:'event'});
      });
    "))
  ),
  footer = div(
    class = "site-footer",
    "NeuroLLPSAI · 云端科研应用原型 · 当前阶段不用于临床诊疗决策"
  ),

  # ====================== 产品首页 ==========================
  tabPanel(
    "首页",
    fluidPage(
      div(
        class = "hero",
        div(class = "hero-kicker", "AI+X · NEUROBLASTOMA PROGNOSIS"),
        h1("NeuroLLPSAI"),
        div(
          class = "hero-sub",
          "面向神经母细胞瘤临床科研与转化研究场景的智能预后评估原型"
        ),
        div(
          class = "hero-mission",
          "将标准化分子数据转化为可解释、可追溯的预后研究信息"
        ),
        div(
          class = "hero-note",
          "数据导入 → 自动质控 → AI推理 → 风险解释 → 报告输出"
        ),
        div(
          class = "hero-actions",
          actionLink("go_assess", "开始智能评估", class = "hero-btn"),
          actionLink("go_report_demo", "查看报告结构", class = "hero-btn-secondary")
        )
      ),

      div(
        class = "metric-row",
        metric_card("493", "GSE49710 主队列"),
        metric_card("7", "核心预后基因"),
        metric_card("0.843", "RSF 内部验证 C-index"),
        metric_card("223", "独立外部验证样本")
      ),

      div(
        class = "card",
        div(class = "section-title", "一个完整的科研转化闭环"),
        div(
          class = "section-sub",
          "NeuroLLPSAI不是单纯展示模型结果，而是把研究成果组织为可操作的分析流程。"
        ),
        div(
          class = "workflow",
          div(class = "flow-card",
              div(class = "flow-num", "STEP 01"),
              div(class = "flow-title", "数据进入"),
              div(class = "flow-text", "上传7基因标准化表达数据，或使用真实验证集演示病例。")),
          div(class = "flow-card",
              div(class = "flow-num", "STEP 02"),
              div(class = "flow-title", "自动质控"),
              div(class = "flow-text", "检查基因完整性、数值格式及训练集表达范围。")),
          div(class = "flow-card",
              div(class = "flow-num", "STEP 03"),
              div(class = "flow-title", "AI推理"),
              div(class = "flow-text", "七基因可解释模型与RSF非线性模型同步评估。")),
          div(class = "flow-card",
              div(class = "flow-num", "STEP 04"),
              div(class = "flow-title", "风险解释"),
              div(class = "flow-text", "输出综合风险方向、风险百分位和1/3/5年预测。")),
          div(class = "flow-card",
              div(class = "flow-num", "STEP 05"),
              div(class = "flow-title", "报告输出"),
              div(class = "flow-text", "形成标准化科研评估报告，便于记录、讨论和后续研究。"))
        )
      ),

      fluidRow(
        column(
          7,
          div(
            class = "card",
            div(class = "section-title", "面向谁使用？"),
            span(class = "audience-tag", "儿童肿瘤临床科研人员"),
            span(class = "audience-tag", "转化医学研究人员"),
            span(class = "audience-tag", "分子检测/生信分析人员"),
            br(), br(),
            div(
              class = "section-sub",
              "当前原型的目标不是让普通患者自行诊断，而是把标准化分子表达结果转化为结构化的研究级风险信息。"
            )
          )
        ),
        column(
          5,
          div(
            class = "card",
            div(class = "section-title", "核心技术"),
            div(class = "feature-text",
                strong("MYCN × LLPS："), "用于定义候选分子特征来源。"),
            br(),
            div(class = "feature-text",
                strong("七基因模型："), "负责可解释的Risk Score与风险分层。"),
            br(),
            div(class = "feature-text",
                strong("RSF："), "在相同七基因输入上学习潜在非线性生存关系。")
          )
        )
      ),

      div(
        class = "card footer-box",
        strong("产品定位："),
        "科研与竞赛展示原型。当前不能替代临床分期、病理、影像、治疗反应或医生综合判断。"
      )
    )
  ),

  # ====================== 智能评估 ==========================
  tabPanel(
    "智能评估",
    fluidPage(
      fluidRow(
        column(
          width = 5,
          div(
            class = "card",
            div(class = "section-title", "01 · 数据输入"),
            div(
              class = "info-strip",
              "建议优先使用“CSV上传”展示标准化数据输入流程；比赛现场也可一键载入真实验证集演示病例。"
            ),
            br(),

            div(
              class = "mode-tabs",
              actionButton("demo_low", "低风险演示"),
              actionButton("demo_high", "高风险演示"),
              downloadButton("download_template", "下载CSV模板")
            ),

            div(
              class = "upload-zone",
              fileInput(
                "csv_file",
                "上传7基因表达CSV",
                accept = c(".csv", "text/csv"),
                buttonLabel = "选择文件",
                placeholder = "未选择文件"
              ),
              textInput(
                "sample_id",
                "样本编号（可匿名）",
                value = "NB-DEMO-001",
                placeholder = "例如 NB-001"
              ),
              div(
                class = "small-note",
                "CSV至少包含两列：Gene 和 Expression。Gene应包含7个目标基因，各出现一次。"
              )
            ),

            br(),
            tags$details(
              class = "advanced-panel",
              tags$summary(
                class = "advanced-summary",
                "高级模式：手动输入7基因表达值"
              ),
              br(),
              uiOutput("manual_inputs"),
              actionButton(
                "use_manual",
                "使用手动输入",
                width = "100%"
              )
            )
          ),

          div(
            class = "card",
            div(class = "section-title", "02 · 数据质控"),
            uiOutput("qc_summary"),
            uiOutput("input_preview"),
            br(),
            actionButton(
              "run_ai",
              "运行AI预后评估",
              class = "btn-primary",
              width = "100%"
            )
          )
        ),

        column(
          width = 7,
          div(
            class = "card",
            div(class = "section-title", "03 · AI推理结果"),
            uiOutput("no_prediction_hint"),
            uiOutput("prediction_area")
          )
        )
      )
    )
  ),

  # ====================== 评估报告 ==========================
  tabPanel(
    "评估报告",
    fluidPage(
      div(
        class = "report-actions no-print",
        actionButton(
          "print_report",
          "打印 / 保存为PDF",
          onclick = "window.print();"
        ),
        downloadButton(
          "download_report",
          "下载HTML评估报告"
        )
      ),
      uiOutput("report_ui")
    )
  ),

  # ====================== 应用与技术 ========================
  tabPanel(
    "应用场景",
    fluidPage(
      div(
        class = "maturity-strip",
        div(class = "maturity-label", "当前阶段"),
        div(
          class = "maturity-text",
          strong("云端可交互科研应用原型。"),
          " 已实现标准化分子数据输入、自动质控、双模型AI推理、风险解释和结构化报告输出；当前用于科研与场景验证，尚未进入临床部署阶段。"
        )
      ),

      div(
        class = "three-grid",
        div(
          class = "feature-card",
          div(class = "feature-title", "科研队列快速分层"),
          div(
            class = "feature-text",
            "将标准化7基因表达数据转化为统一Risk Score、风险方向和RSF预测结果，用于科研队列分层与后续分析。"
          )
        ),
        div(
          class = "feature-card",
          div(class = "feature-title", "转化医学研究"),
          div(
            class = "feature-text",
            "把MYCN背景、LLPS相关特征与生存机器学习结果整合到同一工作流，辅助筛选值得进一步验证的高风险分子特征。"
          )
        ),
        div(
          class = "feature-card",
          div(class = "feature-title", "分子检测结果辅助解释"),
          div(
            class = "feature-text",
            "在研究条件下，将标准化表达结果转化为结构化风险信息，形成可保存、可讨论的科研评估报告。"
          )
        )
      ),

      br(),

      div(
        class = "card",
        div(class = "section-title", "核心技术原理"),
        div(
          class = "section-sub",
          "网页保持应用导向，只展示评委需要理解的核心链路；更详细的方法学、参数和验证过程放在项目文档中。"
        ),

        div(
          class = "tech-principle",

          div(
            class = "tech-stage",
            div(class = "tech-stage-title", "A · 特征发现与可解释建模"),
            div(
              class = "tech-stage-sub",
              "MYCN与LLPS主要作用于上游候选特征定义和模型构建。"
            ),

            div(
              class = "tech-step",
              div(class = "tech-num", "01"),
              div(
                div(class = "tech-step-title", "MYCN相关转录差异"),
                div(
                  class = "tech-step-text",
                  "比较MYCN扩增与非扩增神经母细胞瘤，获得MYCN相关差异表达特征。"
                )
              )
            ),

            div(
              class = "tech-step",
              div(class = "tech-num", "02"),
              div(
                div(class = "tech-step-title", "LLPS候选空间约束"),
                div(
                  class = "tech-step-text",
                  "与DrLLPS相关基因求交，聚焦MYCN相关LLPS分子特征。"
                )
              )
            ),

            div(
              class = "tech-step",
              div(class = "tech-num", "03"),
              div(
                div(class = "tech-step-title", "生存特征筛选"),
                div(
                  class = "tech-step-text",
                  "经Cox与LASSO-Cox等步骤形成七基因核心预后模型。"
                )
              )
            )
          ),

          div(class = "tech-arrow-column", "→"),

          div(
            class = "tech-stage",
            div(class = "tech-stage-title", "B · AI增强与应用推理"),
            div(
              class = "tech-stage-sub",
              "实际使用阶段仅需模型最终确定的7基因标准化表达值。"
            ),

            div(
              class = "tech-step",
              div(class = "tech-num", "04"),
              div(
                div(class = "tech-step-title", "双模型并行推理"),
                div(
                  class = "tech-step-text",
                  "七基因可解释风险模型负责Risk Score与风险分层；RSF在相同输入上学习潜在非线性生存关系。"
                )
              )
            ),

            div(
              class = "tech-step",
              div(class = "tech-num", "05"),
              div(
                div(class = "tech-step-title", "结构化结果输出"),
                div(
                  class = "tech-step-text",
                  "系统将模型结果转化为研究用户更容易使用和记录的风险信息。"
                ),
                div(
                  class = "tech-output-grid",
                  div(class = "tech-output", "综合风险方向"),
                  div(class = "tech-output", "训练集风险百分位"),
                  div(class = "tech-output", "1 / 3 / 5年预测"),
                  div(class = "tech-output", "结构化评估报告")
                )
              )
            )
          )
        ),

        div(
          class = "transparency-note",
          strong("为什么网页没有单独输入MYCN？ "),
          "MYCN并非当前七基因推理阶段的额外输入变量，而是特征发现阶段的生物学锚点；LLPS同样用于候选特征定义。系统实际推理输入为最终确定的7基因标准化表达值。"
        )
      ),

      div(
        class = "evidence-grid",
        div(class = "evidence-card",
            div(class = "evidence-number", "493"),
            div(class = "evidence-label", "GSE49710 主队列")),
        div(class = "evidence-card",
            div(class = "evidence-number", "246 / 247"),
            div(class = "evidence-label", "预定义训练 / 内部验证")),
        div(class = "evidence-card",
            div(class = "evidence-number", "223"),
            div(class = "evidence-label", "E-MTAB-8248 独立外部验证")),
        div(class = "evidence-card",
            div(class = "evidence-number", "RT-qPCR"),
            div(class = "evidence-label", "关键基因表达方向支持"))
      ),

      br(),

      fluidRow(
        column(
          width = 6,
          div(
            class = "card",
            div(class = "section-title", "内部验证：AI推理指标"),
            score_bar("七基因模型 C-index", cindex_values["SevenGene"]),
            score_bar("RSF C-index", cindex_values["RSF"]),
            hr(),
            uiOutput("internal_auc_ui"),
            div(
              class = "small-note",
              "RSF在C-index及1/3/5年AUC上均呈数值提升；当前未进行模型间显著性检验，因此不表述为“显著优于”。"
            )
          )
        ),
        column(
          width = 6,
          div(
            class = "card",
            div(class = "section-title", "RSF变量重要性"),
            uiOutput("vimp_ui"),
            div(
              class = "small-note",
              "Permutation importance用于描述各基因对RSF预测性能的贡献，不等同于Cox回归系数。"
            )
          )
        )
      ),

      fluidRow(
        column(
          width = 6,
          div(
            class = "card",
            div(class = "section-title", "独立外部验证"),
            div(
              class = "section-sub",
              "七基因核心风险模型在E-MTAB-8248（223例）完成独立外部验证。"
            ),
            uiOutput("external_auc_ui"),
            div(
              class = "small-note",
              "RSF增强模块当前主要展示预定义内部验证结果。"
            )
          )
        ),
        column(
          width = 6,
          div(
            class = "card",
            div(class = "section-title", "实验支持"),
            uiOutput("qpcr_ui"),
            br(),
            div(
              class = "small-note",
              "CHD5因表达水平过低/扩增较晚未进行可靠定量；现有实验仅作为表达方向支持，不作直接因果结论。"
            )
          )
        )
      ),

      div(
        class = "card",
        div(class = "section-title", "从科研原型到应用转化"),
        div(
          class = "roadmap",
          div(class = "road-step road-step-current", "当前：云端科研应用原型"),
          div(class = "road-arrow", "→"),
          div(class = "road-step", "多中心回顾性验证"),
          div(class = "road-arrow", "→"),
          div(class = "road-step", "跨平台标准化与校准"),
          div(class = "road-arrow", "→"),
          div(class = "road-step", "前瞻性研究验证"),
          div(class = "road-arrow", "→"),
          div(class = "road-step", "研究级辅助工具")
        ),
        br(),
        div(
          class = "small-note",
          "该路线为项目后续转化规划，不代表当前系统已经具备临床部署资质。"
        )
      )
    )
  )
)

# ------------------------------------------------------------
# 7. Server
# ------------------------------------------------------------
server <- function(input, output, session) {

  # 左上角 NeuroLLPSAI 品牌名：像常规网站一样点击返回首页
  observeEvent(input$brand_home, {
    updateNavbarPage(session, "tabs", selected = "首页")
  }, ignoreInit = TRUE)

  current_values <- reactiveVal(NULL)
  current_source <- reactiveVal("尚未选择数据")
  current_qc <- reactiveVal(NULL)
  prediction <- reactiveVal(NULL)
  last_run_time <- reactiveVal(NULL)

  # 质控动画状态：0=未开始，1~4=依次完成四项检查
  qc_step <- reactiveVal(0)
  qc_running <- reactiveVal(FALSE)
  auto_predict_after_qc <- reactiveVal(FALSE)
  qc_timer <- reactiveTimer(420, session = session)

  observe({
    qc_timer()

    if (!isTRUE(qc_running())) return()

    next_step <- qc_step() + 1
    qc_step(next_step)

    if (next_step >= 4) {
      qc_running(FALSE)

      qc <- current_qc()
      vals <- current_values()

      if (
        isTRUE(auto_predict_after_qc()) &&
        !is.null(qc) &&
        !is.null(vals) &&
        isTRUE(qc$complete_ok) &&
        isTRUE(qc$finite_ok)
      ) {
        prediction(calc_prediction(vals))
        last_run_time(Sys.time())
        showNotification("数据质控完成，AI评估已自动完成。")
      }

      auto_predict_after_qc(FALSE)
    }
  })

  # 首页跳转
  observeEvent(input$go_assess, {
    updateNavbarPage(session, "tabs", selected = "智能评估")
  })

  observeEvent(input$go_report_demo, {
    if (is.null(prediction())) {
      current_values(low_values)
      current_source("真实内部验证集：典型低风险演示病例")
      current_qc(build_qc(low_values, "真实内部验证集演示病例"))
      qc_step(4)
      qc_running(FALSE)
      prediction(calc_prediction(low_values))
      last_run_time(Sys.time())
    }
    updateNavbarPage(session, "tabs", selected = "评估报告")
  })

  # 手动输入框
  output$manual_inputs <- renderUI({
    tagList(
      lapply(
        genes,
        function(g) {
          numericInput(
            paste0("manual_", g),
            paste0(
              g,
              "（训练范围 ",
              round(gene_min[g], 2),
              "–",
              round(gene_max[g], 2),
              "）"
            ),
            value = round(low_values[g], 5),
            step = 0.001
          )
        }
      )
    )
  })

  # 下载CSV模板
  output$download_template <- downloadHandler(
    filename = function() "NeuroLLPSAI_7gene_template.csv",
    content = function(file) {
      df <- data.frame(
        Gene = genes,
        Expression = round(low_values[genes], 5),
        stringsAsFactors = FALSE
      )
      write.csv(df, file, row.names = FALSE)
    }
  )

  set_case <- function(values, source_name, sample_id = NULL, auto_predict = FALSE) {
    values <- as.numeric(values)
    names(values) <- genes

    current_values(values)
    current_source(source_name)
    current_qc(build_qc(values, source_name))

    # 每次换数据源都从“待检测”重新开始，逐项完成质控
    qc_step(0)
    qc_running(TRUE)
    auto_predict_after_qc(auto_predict)
    prediction(NULL)

    if (!is.null(sample_id) && nzchar(sample_id)) {
      updateTextInput(session, "sample_id", value = sample_id)
    }
  }

  # 演示病例：自动完成评估
  observeEvent(input$demo_low, {
    set_case(
      low_values,
      "真实内部验证集：典型低风险演示病例",
      sample_id = "NB-DEMO-LOW",
      auto_predict = TRUE
    )
    showNotification("已载入低风险演示病例，正在进行自动质控。")
  })

  observeEvent(input$demo_high, {
    set_case(
      high_values,
      "真实内部验证集：典型高风险演示病例",
      sample_id = "NB-DEMO-HIGH",
      auto_predict = TRUE
    )
    showNotification("已载入高风险演示病例，正在进行自动质控。")
  })

  # CSV上传
  observeEvent(input$csv_file, {
    req(input$csv_file)

    df <- tryCatch(
      read.csv(
        input$csv_file$datapath,
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),
      error = function(e) NULL
    )

    if (is.null(df)) {
      current_values(NULL)
      qc_step(0)
      qc_running(FALSE)
      auto_predict_after_qc(FALSE)
      current_qc( list(
        source = paste0("CSV：", input$csv_file$name),
        complete_ok = FALSE,
        finite_ok = FALSE,
        range_ok = FALSE,
        outside = character(0),
        file_error = TRUE,
        message = "CSV读取失败。"
      ))
      prediction(NULL)
      return()
    }

    cn <- tolower(trimws(colnames(df)))
    gene_col <- which(cn == "gene")
    expr_col <- which(cn == "expression")

    if (length(gene_col) != 1 || length(expr_col) != 1) {
      current_values(NULL)
      qc_step(0)
      qc_running(FALSE)
      auto_predict_after_qc(FALSE)
      current_qc(list(
        source = paste0("CSV：", input$csv_file$name),
        complete_ok = FALSE,
        finite_ok = FALSE,
        range_ok = FALSE,
        outside = character(0),
        file_error = TRUE,
        message = "CSV必须包含 Gene 和 Expression 两列。"
      ))
      prediction(NULL)
      return()
    }

    g <- toupper(trimws(as.character(df[[gene_col]])))
    e <- suppressWarnings(as.numeric(df[[expr_col]]))

    duplicate_genes <- unique(g[duplicated(g)])
    required_present <- genes %in% g
    complete_ok <- all(required_present) && length(duplicate_genes) == 0

    if (!complete_ok) {
      msg <- c()
      missing <- genes[!required_present]
      if (length(missing) > 0) {
        msg <- c(msg, paste0("缺少：", paste(missing, collapse = "、")))
      }
      if (length(duplicate_genes) > 0) {
        msg <- c(msg, paste0("重复：", paste(duplicate_genes, collapse = "、")))
      }

      current_values(NULL)
      qc_step(0)
      qc_running(FALSE)
      auto_predict_after_qc(FALSE)
      current_qc(list(
        source = paste0("CSV：", input$csv_file$name),
        complete_ok = FALSE,
        finite_ok = FALSE,
        range_ok = FALSE,
        outside = character(0),
        file_error = TRUE,
        message = paste(msg, collapse = "；")
      ))
      prediction(NULL)
      return()
    }

    values <- setNames(
      sapply(
        genes,
        function(x) e[match(x, g)]
      ),
      genes
    )

    finite_ok <- all(is.finite(values))
    outside <- genes[
      is.finite(values) &
      (values < gene_min[genes] | values > gene_max[genes])
    ]

    qc <- list(
      source = paste0("CSV：", input$csv_file$name),
      complete_ok = complete_ok,
      finite_ok = finite_ok,
      range_ok = finite_ok && length(outside) == 0,
      outside = outside,
      file_error = FALSE,
      message = NULL
    )

    current_values(values)
    current_source(qc$source)
    current_qc(qc)
    qc_step(0)
    qc_running(TRUE)
    auto_predict_after_qc(FALSE)
    prediction(NULL)

    showNotification("CSV读取完成，正在逐项进行数据质控。")
  })

  # 手动输入
  observeEvent(input$use_manual, {
    values <- sapply(
      genes,
      function(g) input[[paste0("manual_", g)]]
    )
    names(values) <- genes

    if (!all(is.finite(values))) {
      showNotification("请确保7个基因均填写有效数字。", type = "error")
      return()
    }

    set_case(
      values,
      "高级模式：手动输入",
      auto_predict = FALSE
    )

    showNotification("已采用手动输入，正在逐项进行数据质控。")
  })

  # 运行AI
  observeEvent(input$run_ai, {
    qc <- current_qc()
    vals <- current_values()

    if (isTRUE(qc_running()) || qc_step() < 4) {
      showNotification("数据质控尚未完成，请稍候。", type = "warning")
      return()
    }

    if (is.null(qc) || !isTRUE(qc$complete_ok) || !isTRUE(qc$finite_ok)) {
      showNotification("当前数据未通过必要质控，无法运行AI评估。", type = "error")
      return()
    }

    prediction(calc_prediction(vals))
    last_run_time(Sys.time())
    showNotification("AI评估完成，可查看结果或生成报告。")
  })

  # 数据质控
  output$qc_summary <- renderUI({
    qc <- current_qc()
    step <- qc_step()

    pending_badge <- function(label) {
      div(
        class = "qc-badge qc-pending",
        div(class = "qc-symbol", "○"),
        div(class = "qc-text", label)
      )
    }

    staged_badge <- function(index, ok, good, bad, warn = FALSE) {
      if (step < index) {
        return(pending_badge(paste0("待检测：", good)))
      }
      status_badge(ok, good, bad, warn = warn)
    }

    # 初次进入页面：不预设任何患者数据，也不提前打勾
    if (is.null(qc)) {
      return(
        tagList(
          div(class = "small-note", "当前数据源：尚未选择"),
          div(
            class = "qc-grid",
            pending_badge("待检测：数据格式"),
            pending_badge("待检测：7/7基因"),
            pending_badge("待检测：表达数值"),
            pending_badge("待检测：训练范围")
          ),
          div(
            class = "small-note",
            "选择演示病例、上传CSV或使用手动输入后，系统会依次执行四项检查。"
          )
        )
      )
    }

    if (isTRUE(qc$file_error)) {
      return(div(
        class = "warn-strip",
        strong("数据读取未通过："),
        qc$message
      ))
    }

    div(
      tagList(
        div(class = "small-note", paste0("当前数据源：", qc$source)),
        div(
          class = "qc-grid",
          staged_badge(
            1,
            TRUE,
            "数据格式可读取",
            "格式错误"
          ),
          staged_badge(
            2,
            qc$complete_ok,
            "7/7基因完整",
            "基因不完整"
          ),
          staged_badge(
            3,
            qc$finite_ok,
            "表达值均为数值",
            "存在无效数值"
          ),
          staged_badge(
            4,
            qc$range_ok,
            "均在训练范围",
            if (length(qc$outside) > 0) {
              paste0("超范围：", paste(qc$outside, collapse = "、"))
            } else {
              "范围需检查"
            },
            warn = !qc$range_ok && qc$finite_ok
          )
        ),
        if (isTRUE(qc_running())) {
          div(class = "small-note", "正在逐项检查数据，请稍候……")
        } else {
          div(class = "small-note", "数据质控完成。通过必要检查后即可运行AI预后评估。")
        }
      )
    )
  })

  output$input_preview <- renderUI({
    vals <- current_values()
    if (is.null(vals)) return(NULL)

    rows <- lapply(
      genes,
      function(g) {
        tags$tr(
          tags$td(g),
          tags$td(sprintf("%.5f", vals[g])),
          tags$td(
            paste0(
              round(gene_min[g], 2),
              " – ",
              round(gene_max[g], 2)
            )
          )
        )
      }
    )

    tagList(
      tags$table(
        class = "data-table",
        tags$thead(
          tags$tr(
            tags$th("Gene"),
            tags$th("Expression"),
            tags$th("训练范围")
          )
        ),
        tags$tbody(rows)
      )
    )
  })

  # 预测结果
  output$no_prediction_hint <- renderUI({
    if (!is.null(prediction())) return(NULL)

    div(
      class = "info-strip",
      strong("等待AI推理："),
      "当前数据已进入系统。通过必要质控后点击“运行AI预后评估”。"
    )
  })

  output$prediction_area <- renderUI({
    x <- prediction()
    if (is.null(x)) return(NULL)

    overall <- overall_text(x)

    tagList(
      div(
        class = "result-hero",
        div(class = "result-eyebrow", "双模型综合结果"),
        div(class = "result-main", paste0("综合判断：", overall$level)),
        div(class = "result-sub", overall$brief)
      ),

      div(
        class = "model-grid",
        div(
          class = "model-card",
          div(class = "model-name", "七基因可解释风险模型"),
          div(
            class = "model-judgement",
            if (x$cox_high) "高风险" else "低风险"
          ),
          div(
            class = "small-note",
            paste0("Risk Score = ", sprintf("%.3f", x$cox_score))
          ),
          risk_gauge(x$cox_pct)
        ),
        div(
          class = "model-card",
          div(class = "model-name", "随机生存森林（RSF）增强模型"),
          div(
            class = "model-judgement",
            if (x$rsf_high) "高风险倾向" else "低风险倾向"
          ),
          div(
            class = "small-note",
            paste0("RSF 风险评分 = ", sprintf("%.3f", x$rsf_score))
          ),
          risk_gauge(x$rsf_pct)
        )
      ),

      div(
        class = "survival-summary",
        div(class = "surv-mini",
            div(class = "surv-year", "1年预测生存概率"),
            div(class = "surv-value", sprintf("%.1f%%", 100 * x$surv_1y))),
        div(class = "surv-mini",
            div(class = "surv-year", "3年预测生存概率"),
            div(class = "surv-value", sprintf("%.1f%%", 100 * x$surv_3y))),
        div(class = "surv-mini",
            div(class = "surv-year", "5年预测生存概率"),
            div(class = "surv-value", sprintf("%.1f%%", 100 * x$surv_5y)))
      ),

      plotOutput("assessment_survival_plot", height = "290px"),

      div(
        class = "good-strip",
        strong("下一步："),
        actionLink("go_report", "生成结构化评估报告")
      )
    )
  })

  output$assessment_survival_plot <- renderPlot({
    x <- prediction()
    req(x)

    years <- x$pred_times / 365

    par(mar = c(4.5, 4.5, 2.5, 1.5))

    plot(
      years,
      x$surv_curve,
      type = "s",
      lwd = 3,
      col = "#1e638f",
      xlab = "随访时间（年）",
      ylab = "预测生存概率",
      ylim = c(0, 1),
      xlim = c(0, min(8, max(years, na.rm = TRUE))),
      main = "RSF个体预测生存曲线",
      bty = "l",
      las = 1
    )

    abline(v = c(1, 3, 5), lty = 3, col = "#9ba7b2")
    points(
      c(1, 3, 5),
      c(x$surv_1y, x$surv_3y, x$surv_5y),
      pch = 19,
      cex = 1.15,
      col = "#102f4f"
    )
  })

  observeEvent(input$go_report, {
    updateNavbarPage(session, "tabs", selected = "评估报告")
  })

  # 报告页面
  output$report_ui <- renderUI({
    x <- prediction()

    if (is.null(x)) {
      return(
        div(
          class = "card",
          div(class = "section-title", "尚未生成评估报告"),
          div(
            class = "section-sub",
            "请先在“智能评估”页面载入或上传数据，并完成AI预后评估。"
          ),
          actionLink("back_assess", "前往智能评估")
        )
      )
    }

    overall <- overall_text(x)
    sid <- if (!is.null(input$sample_id) && nzchar(input$sample_id)) {
      input$sample_id
    } else {
      "未命名样本"
    }

    time_text <- if (!is.null(last_run_time())) {
      format(last_run_time(), "%Y-%m-%d %H:%M")
    } else {
      format(Sys.time(), "%Y-%m-%d %H:%M")
    }

    tagList(
      div(
        class = "report-shell",

        div(
          class = "report-header",
          div(
            div(class = "report-title", "NeuroLLPSAI 智能预后评估报告"),
            div(class = "small-note", "科研与竞赛展示原型")
          ),
          div(
            class = "report-meta",
            div(paste0("样本编号：", sid)),
            div(paste0("数据来源：", current_source())),
            div(paste0("生成时间：", time_text))
          )
        ),

        div(
          class = "report-section",
          div(class = "report-section-title", "1. 综合风险结论"),
          div(class = "result-main", overall$level),
          div(class = "result-sub", overall$brief)
        ),

        div(
          class = "report-section",
          div(class = "report-section-title", "2. 双模型结果"),
          div(
            class = "model-grid",
            div(
              class = "model-card",
              div(class = "model-name", "七基因可解释风险模型"),
              div(
                class = "model-judgement",
                if (x$cox_high) "高风险" else "低风险"
              ),
              div(class = "small-note",
                  paste0("Risk Score = ", sprintf("%.3f", x$cox_score))),
              risk_gauge(x$cox_pct)
            ),
            div(
              class = "model-card",
              div(class = "model-name", "随机生存森林（RSF）增强模型"),
              div(
                class = "model-judgement",
                if (x$rsf_high) "高风险倾向" else "低风险倾向"
              ),
              div(class = "small-note",
                  paste0("RSF 风险评分 = ", sprintf("%.3f", x$rsf_score))),
              risk_gauge(x$rsf_pct)
            )
          )
        ),

        div(
          class = "report-section",
          div(class = "report-section-title", "3. 预测生存概率"),
          div(
            class = "survival-summary",
            div(class = "surv-mini",
                div(class = "surv-year", "1年"),
                div(class = "surv-value", sprintf("%.1f%%", 100 * x$surv_1y))),
            div(class = "surv-mini",
                div(class = "surv-year", "3年"),
                div(class = "surv-value", sprintf("%.1f%%", 100 * x$surv_3y))),
            div(class = "surv-mini",
                div(class = "surv-year", "5年"),
                div(class = "surv-value", sprintf("%.1f%%", 100 * x$surv_5y)))
          )
        ),

        div(
          class = "report-section",
          div(class = "report-section-title", "4. 研究解释"),
          if (overall$class == "high") {
            div(
              class = "warn-strip",
              "该样本在当前研究模型中呈较高预后风险方向。建议在科研讨论中结合年龄、分期、MYCN状态、治疗反应等完整临床信息综合解释。"
            )
          } else if (overall$class == "low") {
            div(
              class = "good-strip",
              "该样本在当前研究模型中呈较低预后风险方向。该结果仅代表模型推理，应与年龄、分期、MYCN状态、治疗反应等完整临床信息共同解释。"
            )
          } else {
            div(
              class = "warn-strip",
              "两种模型对该样本的风险方向判断不完全一致，提示应特别谨慎解释，并结合完整临床资料进行研究判断。"
            )
          }
        ),

        div(
          class = "report-section",
          div(class = "report-section-title", "5. 使用边界"),
          div(
            class = "small-note",
            "本报告由科研原型自动生成，不构成临床诊断、治疗决策或个体化医疗建议。当前模型主要基于公开回顾性队列，实际转化仍需多中心、跨平台和前瞻性验证。"
          )
        )
      )
    )
  })

  observeEvent(input$back_assess, {
    updateNavbarPage(session, "tabs", selected = "智能评估")
  })

  # 下载HTML报告
  output$download_report <- downloadHandler(
    filename = function() {
      sid <- if (!is.null(input$sample_id) && nzchar(input$sample_id)) {
        gsub("[^A-Za-z0-9_-]", "_", input$sample_id)
      } else {
        "sample"
      }
      paste0("NeuroLLPSAI_", sid, "_report.html")
    },
    content = function(file) {
      x <- prediction()

      if (is.null(x)) {
        writeLines("<html><body><p>尚未完成AI评估。</p></body></html>", file)
        return()
      }

      overall <- overall_text(x)
      sid <- if (!is.null(input$sample_id) && nzchar(input$sample_id)) {
        input$sample_id
      } else {
        "未命名样本"
      }

      gene_rows <- paste0(
        "<tr><td>", genes, "</td><td>",
        sprintf("%.5f", x$values[genes]),
        "</td></tr>",
        collapse = "\n"
      )

      html <- paste0(
'<!doctype html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<title>NeuroLLPSAI Report</title>
<style>
body{font-family:Arial,"Microsoft YaHei",sans-serif;margin:40px;color:#172033;line-height:1.7;}
h1{color:#164e72;} .box{border:1px solid #dfe7ee;border-radius:12px;padding:16px;margin:14px 0;}
.big{font-size:28px;font-weight:800;} table{border-collapse:collapse;width:100%;}
th,td{padding:8px;border-bottom:1px solid #e7edf3;text-align:left;} .note{color:#687786;font-size:12px;}
</style>
</head>
<body>
<h1>NeuroLLPSAI 智能预后评估报告</h1>
<p>样本编号：', safe_html(sid), '</p>
<p>数据来源：', safe_html(current_source()), '</p>
<div class="box"><div class="big">', safe_html(overall$level), '</div><p>',
safe_html(overall$brief), '</p></div>
<div class="box">
<h3>双模型结果</h3>
<p>七基因模型：', if (x$cox_high) "高风险" else "低风险",
'；Risk Score = ', sprintf("%.3f", x$cox_score),
'；训练集风险位置第 ', round(x$cox_pct), ' 百分位。</p>
<p>RSF模型：', if (x$rsf_high) "高风险倾向" else "低风险倾向",
'；RSF风险评分 = ', sprintf("%.3f", x$rsf_score),
'；训练集风险位置第 ', round(x$rsf_pct), ' 百分位。</p>
</div>
<div class="box">
<h3>预测生存概率</h3>
<p>1年：', sprintf("%.1f%%", 100*x$surv_1y),
'　3年：', sprintf("%.1f%%", 100*x$surv_3y),
'　5年：', sprintf("%.1f%%", 100*x$surv_5y), '</p>
</div>
<div class="box">
<h3>输入的7基因表达值</h3>
<table><tr><th>Gene</th><th>Expression</th></tr>', gene_rows, '</table>
</div>
<p class="note">本报告由科研与竞赛展示原型自动生成，不用于临床诊断、治疗决策或个体化医疗建议。</p>
</body></html>'
      )

      writeLines(html, file, useBytes = TRUE)
    }
  )

  # 应用与技术：AUC比较
  output$internal_auc_ui <- renderUI({
    tagList(
      lapply(
        seq_len(nrow(internal_perf)),
        function(i) {
          row <- internal_perf[i, ]
          div(
            class = "pair-row",
            div(class = "pair-time", row$Time),
            div(
              class = "pair-values",
              div(
                class = "pair-chip",
                span("七基因模型"),
                strong(sprintf("%.3f", row$SevenGene))
              ),
              div(
                class = "pair-chip",
                span("RSF"),
                strong(sprintf("%.3f", row$RSF))
              )
            )
          )
        }
      )
    )
  })

  output$vimp_ui <- renderUI({
    vmax <- max(vimp$Importance, na.rm = TRUE)

    tagList(
      lapply(
        seq_len(nrow(vimp)),
        function(i) {
          pct <- 100 * vimp$Importance[i] / vmax
          div(
            class = "vimp-row",
            div(
              class = "vimp-head",
              span(vimp$Gene[i]),
              span(sprintf("%.3f", vimp$Importance[i]))
            ),
            div(
              class = "vimp-track",
              div(
                class = "vimp-fill",
                style = paste0("width:", pct, "%;")
              )
            )
          )
        }
      )
    )
  })

  output$external_auc_ui <- renderUI({
    tagList(
      lapply(
        names(external_auc),
        function(nm) {
          score_bar(
            paste0(nm, " AUC"),
            external_auc[[nm]]
          )
        }
      )
    )
  })

  output$qpcr_ui <- renderUI({
    tagList(
      lapply(
        names(qpcr_result),
        function(g) {
          div(
            class = "pair-row",
            div(class = "pair-time", g),
            div(class = "pair-values",
                div(class = "pair-chip", qpcr_result[[g]]))
          )
        }
      )
    )
  })
}

shinyApp(ui = ui, server = server)
