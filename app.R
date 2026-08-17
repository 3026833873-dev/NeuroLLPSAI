# ============================================================
# NeuroLLPSAI Web Test
# 云端发布测试版
# 研究与竞赛展示用途，不用于临床诊疗决策
# ============================================================

library(shiny)
library(randomForestSRC)

# 云端使用相对路径，不依赖任何电脑的C盘
paper_model_file <- file.path("model", "seven_gene_model.rds")
rsf_model_file <- file.path("model", "rsf_model.rds")

if (!file.exists(paper_model_file)) {
  stop("找不到 model/seven_gene_model.rds")
}
if (!file.exists(rsf_model_file)) {
  stop("找不到 model/rsf_model.rds")
}

paper <- readRDS(paper_model_file)
rsf_model <- readRDS(rsf_model_file)

genes <- paper$genes
model_coef <- paper$model_coef[genes]
risk_cutoff <- paper$risk_cutoff
expr_train <- paper$expr_train
expr_valid <- paper$expr_valid
train_risk_df <- paper$train_risk_df
valid_risk_df <- paper$valid_risk_df

gene_min <- apply(expr_train, 1, min, na.rm = TRUE)
gene_max <- apply(expr_train, 1, max, na.rm = TRUE)

rsf_train_score <- as.numeric(rsf_model$predicted)
rsf_cutoff <- median(rsf_train_score, na.rm = TRUE)

# 从247例预定义内部验证集中选择稳健演示病例
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

cox_sd <- sd(cox_valid_score, na.rm = TRUE)
rsf_sd <- sd(rsf_valid_score, na.rm = TRUE)
if (!is.finite(cox_sd) || cox_sd <= 0) cox_sd <- 1
if (!is.finite(rsf_sd) || rsf_sd <= 0) rsf_sd <- 1

low_candidates <- which(!cox_valid_high & !rsf_valid_high & within_training_range)
high_candidates <- which(cox_valid_high & rsf_valid_high & within_training_range)

if (length(low_candidates) == 0 || length(high_candidates) == 0) {
  stop("未找到适合的验证集演示病例。")
}

low_strength <- pmin(
  (risk_cutoff - cox_valid_score[low_candidates]) / cox_sd,
  (rsf_cutoff - rsf_valid_score[low_candidates]) / rsf_sd
)

high_strength <- pmin(
  (cox_valid_score[high_candidates] - risk_cutoff) / cox_sd,
  (rsf_valid_score[high_candidates] - rsf_cutoff) / rsf_sd
)

low_idx <- low_candidates[which.max(low_strength)]
high_idx <- high_candidates[which.max(high_strength)]

low_values <- expr_valid[, valid_samples[low_idx]]
high_values <- expr_valid[, valid_samples[high_idx]]
default_values <- low_values

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body {background:#f5f7fa;font-family:'Microsoft YaHei',Arial,sans-serif;}
      .title-box,.panel-box,.result-card {
        background:white;border-radius:12px;box-shadow:0 2px 10px rgba(0,0,0,.06);
      }
      .title-box {padding:20px 24px;margin-bottom:18px;}
      .panel-box {padding:18px;margin-bottom:16px;}
      .result-card {padding:16px;margin-bottom:12px;border-left:5px solid #4f6bed;}
      .big-number {font-size:28px;font-weight:700;margin:4px 0;}
      .small-note {color:#666;font-size:12px;line-height:1.6;}
      .demo-note {background:#eef4ff;padding:10px 12px;border-radius:8px;margin-bottom:12px;}
      .disclaimer {background:#fff8e6;padding:12px 16px;border-radius:8px;margin-top:14px;}
    "))
  ),

  div(
    class="title-box",
    h2("NeuroLLPSAI 神经母细胞瘤智能预后预测演示系统"),
    p("融合 MYCN 与液-液相分离相关七基因特征的可解释模型 + 随机生存森林（RSF）"),
    tags$span(class="small-note","云端测试版｜研究与竞赛展示原型｜不用于临床诊疗决策")
  ),

  fluidRow(
    column(
      4,
      div(
        class="panel-box",
        h4("① 输入7个基因表达值"),
        div(class="demo-note", strong("演示建议："), "直接载入验证集示例。"),
        actionButton("load_low","载入验证集低风险示例",width="100%"),
        br(),br(),
        actionButton("load_high","载入验证集高风险示例",width="100%"),
        br(),br(),
        lapply(genes,function(g){
          numericInput(
            paste0("gene_",g),
            paste0(g,"（训练范围 ",round(gene_min[g],2),"–",round(gene_max[g],2),"）"),
            value=round(as.numeric(default_values[g]),5),
            step=0.001
          )
        }),
        actionButton("predict_btn","开始预测",class="btn-primary",width="100%")
      )
    ),

    column(
      8,
      div(
        class="panel-box",
        h4("② 预测结果"),
        uiOutput("warning_box"),
        uiOutput("cox_result"),
        uiOutput("rsf_result"),
        uiOutput("survival_result"),
        uiOutput("agreement_result")
      ),
      div(
        class="panel-box",
        h4("③ 预定义内部验证集中的模型表现"),
        tableOutput("performance_table"),
        tags$span(
          class="small-note",
          "RSF在本验证集中呈现数值上的性能提升；当前未进行模型间显著性检验，因此不表述为“显著优于”。"
        )
      )
    )
  ),

  div(
    class="disclaimer",
    strong("重要说明："),
    "本系统为科研与竞赛展示原型。输入值应采用与训练数据一致的表达数据处理方式；预测结果不能替代临床综合判断。"
  )
)

server <- function(input, output, session) {

  observeEvent(input$load_low, {
    for (g in genes) {
      updateNumericInput(session,paste0("gene_",g),value=round(as.numeric(low_values[g]),5))
    }
  })

  observeEvent(input$load_high, {
    for (g in genes) {
      updateNumericInput(session,paste0("gene_",g),value=round(as.numeric(high_values[g]),5))
    }
  })

  prediction <- eventReactive(input$predict_btn, {
    values <- sapply(genes,function(g) input[[paste0("gene_",g)]])
    names(values) <- genes

    validate(need(all(is.finite(values)),"请确保7个基因都填写了有效数字。"))

    outside <- genes[values < gene_min[genes] | values > gene_max[genes]]

    cox_score <- sum(values * model_coef[genes])
    cox_high <- cox_score > risk_cutoff

    newdata <- as.data.frame(as.list(values))
    newdata <- newdata[,genes,drop=FALSE]
    rsf_pred <- predict(rsf_model,newdata=newdata)
    rsf_score <- as.numeric(rsf_pred$predicted)[1]
    rsf_high <- rsf_score > rsf_cutoff

    pred_times <- rsf_pred$time.interest
    surv_curve <- rsf_pred$survival
    if (!is.null(dim(surv_curve))) surv_curve <- surv_curve[1,]
    surv_curve <- as.numeric(surv_curve)

    get_surv <- function(day){
      idx <- which.min(abs(pred_times-day))
      max(min(surv_curve[idx],1),0)
    }

    list(
      outside=outside,
      cox_score=cox_score,
      cox_high=cox_high,
      rsf_score=rsf_score,
      rsf_high=rsf_high,
      surv_1y=get_surv(365),
      surv_3y=get_surv(1095),
      surv_5y=get_surv(1825)
    )
  })

  output$warning_box <- renderUI({
    req(prediction())
    x <- prediction()
    if (length(x$outside)==0) return(NULL)
    div(class="disclaimer",strong("输入范围提示："),
        paste0("以下基因超出训练集观察范围：",paste(x$outside,collapse="、"),"。"))
  })

  output$cox_result <- renderUI({
    req(prediction()); x <- prediction()
    div(class="result-card",
        h4("七基因可解释风险模型"),
        div(class="big-number",if(x$cox_high)"高风险" else "低风险"),
        p(paste0("Risk Score = ",sprintf("%.3f",x$cox_score))),
        tags$span(class="small-note",
                  paste0("训练集风险分界值：",sprintf("%.3f",risk_cutoff),"。")))
  })

  output$rsf_result <- renderUI({
    req(prediction()); x <- prediction()
    div(class="result-card",
        h4("随机生存森林（RSF）增强模型"),
        div(class="big-number",if(x$rsf_high)"高风险倾向" else "低风险倾向"),
        p(paste0("RSF 风险评分 = ",sprintf("%.3f",x$rsf_score))))
  })

  output$survival_result <- renderUI({
    req(prediction()); x <- prediction()
    div(class="result-card",
        h4("RSF预测生存概率（研究演示）"),
        fluidRow(
          column(4,strong("1年"),div(class="big-number",sprintf("%.1f%%",100*x$surv_1y))),
          column(4,strong("3年"),div(class="big-number",sprintf("%.1f%%",100*x$surv_3y))),
          column(4,strong("5年"),div(class="big-number",sprintf("%.1f%%",100*x$surv_5y)))
        ))
  })

  output$agreement_result <- renderUI({
    req(prediction()); x <- prediction()
    agree <- identical(x$cox_high,x$rsf_high)
    txt <- if (agree) "两种模型提示的风险方向一致。" else
      "两种模型的风险方向不完全一致，应结合完整临床信息谨慎解释。"
    div(class="result-card",h4("双模型结果一致性"),p(txt))
  })

  output$performance_table <- renderTable(
    data.frame(
      模型=c("七基因可解释模型","随机生存森林（RSF）"),
      `C-index`=c(0.832,0.843),
      `1年AUC`=c(0.758,0.798),
      `3年AUC`=c(0.871,0.880),
      `5年AUC`=c(0.882,0.888),
      check.names=FALSE
    ),
    digits=3,striped=TRUE,bordered=TRUE
  )
}

shinyApp(ui=ui,server=server)
