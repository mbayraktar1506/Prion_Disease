# Gerekli Paketlerin Yüklenmesi
if (!requireNamespace("BiocManager",  quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("GEOquery",     quietly = TRUE)) BiocManager::install("GEOquery")
if (!requireNamespace("DESeq2",       quietly = TRUE)) BiocManager::install("DESeq2")
if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) BiocManager::install("org.Mm.eg.db")
if (!requireNamespace("ggplot2",      quietly = TRUE)) install.packages("ggplot2")

library(GEOquery)
library(DESeq2)
library(org.Mm.eg.db)
library(ggplot2)





# Hasta/Örnek Bilgilerini (Metadata) Çekme
cat("Örnek bilgileri GEO'dan indiriliyor...\n")
gse         <- getGEO("GSE198063", AnnotGPL = TRUE)[[1]]
ornek_bilgi <- pData(gse)



# Ham (Raw) Count Verisini İndirme
cat("Ham (raw count) veri dosyası indiriliyor...\n")
link_raw <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198063/suppl/GSE198063_raw_count_table.csv.gz"
download.file(link_raw, "GSE198063_raw_count.csv.gz", mode = "wb")
raw_veri <- read.csv("GSE198063_raw_count.csv.gz", row.names = 1, check.names = FALSE)




#Eşleştirme düzeltmesi
kodlar <- gsub(".", "-", colnames(raw_veri), fixed = TRUE)
sira   <- match(kodlar, ornek_bilgi$description)
cat("Eşleşmeyen örnek sayısı:", sum(is.na(sira)), "(0 olmalı)\n")

if (sum(is.na(sira)) > 0) {
  stop("Bazı örnekler GEO metadata'sıyla eşleşmedi! colnames(raw_veri) ve ",
       "ornek_bilgi$description alanlarını kontrol edin.")
}

meta <- ornek_bilgi[sira, ]
rownames(meta) <- colnames(raw_veri)   # artık meta ile raw_veri aynı sırada




#vGluT2 + Cerebrum + RiboTag-IP Alt Kümesinin Seçilmesi
genotip  <- meta[["genotype:ch1"]]
hucre    <- meta[["cell type:ch1"]]
bolge    <- meta[["source_name_ch1"]]
protokol <- meta[["protocol:ch1"]]

sec <- grepl("RiboTag",       protokol, ignore.case = TRUE) &
  grepl("glut|slc17a6",  hucre,    ignore.case = TRUE) &
  grepl("Cerebrum",      bolge,    ignore.case = TRUE) &
  !grepl("BE268", colnames(raw_veri), ignore.case = TRUE)


cat("Seçilen (vGluT2 + Cerebrum + IP) örnek sayısı:", sum(sec), "\n")

raw_veri_alt <- raw_veri[, sec]
genotip_alt  <- genotip[sec]



#Hastaları Genotipe Göre Gruplama (CJD / FFI / Control)
grup <- rep(NA, length(genotip_alt))
grup[grepl("Control|WT", genotip_alt, ignore.case = TRUE)] <- "Control"
grup[grepl("CJD",        genotip_alt, ignore.case = TRUE)] <- "CJD"
grup[grepl("FFI",        genotip_alt, ignore.case = TRUE)] <- "FFI"
cat("Grup dağılımı:\n"); print(table(grup, useNA = "always"))



# Grubu belirlenemeyen (NA) örnek kalırsa at:
gecerli      <- !is.na(grup)
raw_veri_alt <- raw_veri_alt[, gecerli, drop = FALSE]
grup         <- grup[gecerli]





#DESeq2 Analiz Motorunun Kurulması
cat("DESeq2 analizi başlatılıyor...\n")

coldata <- data.frame(condition = factor(grup, levels = c("Control", "CJD", "FFI")))
rownames(coldata) <- colnames(raw_veri_alt)   # sıra zaten garanti, isimler de eşleşiyor



# Count matrisi tam sayı olmalı (raw veri zaten ham sayım, güvenlik için round):
raw_veri_alt <- round(as.matrix(raw_veri_alt))
storage.mode(raw_veri_alt) <- "integer"

dds <- DESeqDataSetFromMatrix(countData = raw_veri_alt,
                              colData   = coldata,
                              design    = ~ condition)


#Ön Filtreleme (makaledeki kriter:ortalama sayım > 10)
keep <- rowMeans(counts(dds)) > 10
dds  <- dds[keep, ]
cat("Ön filtreden sonra kalan gen sayısı:", nrow(dds), "\n")



# Makaledeki DESeq2 normalizasyonu (size factor + dispersion + Wald testi):
dds <- DESeq(dds)




#Normalizasyon Öncesi/Sonrası Karşılaştırma + PCA
cat("Normalizasyon grafikleri ve PCA çiziliyor...\n")

vsd <- vst(dds, blind = FALSE)

png("GSE198063_vGluT2_before_after_normalizasyon.png", width = 1000, height = 500)
par(mfrow = c(1, 2))
boxplot(log2(counts(dds, normalized = FALSE) + 1),
        main = "Before Normalization (Raw, log2)", las = 2, outline = FALSE,
        col = c("Control" = "steelblue", "CJD" = "tomato", "FFI" = "forestgreen")[as.character(coldata$condition)])
boxplot(assay(vsd),
        main = "After DESeq2 Normalization (VST)", las = 2, outline = FALSE,
        col = c("Control" = "steelblue", "CJD" = "tomato", "FFI" = "forestgreen")[as.character(coldata$condition)])
par(mfrow = c(1, 1))
dev.off()

#PCA: NORMALİZASYON SONRASI (DESeq2 VST üzerinden) ---
pca_grafik <- plotPCA(vsd, intgroup = "condition") +
  geom_text(aes(label = name), vjust = -0.7, size = 3, show.legend = FALSE) +
  ggtitle("PCA - Normalizasyon SONRASI (DESeq2 VST)") +
  theme_minimal()
ggsave("GSE198063_vGluT2_PCA_normalizasyon_sonrasi.png", pca_grafik, width = 7, height = 6)
print(pca_grafik)




#Gen Sembolü Eşleme Fonksiyonu (Ensembl -> Symbol, fare)
sembol_getir <- function(ids) {
  if (length(ids) == 0) return(character(0))
  
  # Bazı GEO/Ensembl dosyalarında ID'lerin sonunda versiyon eki varsa temizle
  ids_temiz <- sub("\\.[0-9]+$", "", ids)
  
  mapIds(org.Mm.eg.db, keys = ids_temiz, keytype = "ENSEMBL",
         column = "SYMBOL", multiVals = "first")
}



#DEG Çıkarma Fonksiyonu
sonuc_cikar_ve_kaydet <- function(grup_A, grup_B, p_siniri = 0.05, lfc_siniri = NULL,
                                  kriter = c("padj", "pvalue")) {
  kriter <- match.arg(kriter)
  cat("\n--->", grup_A, "vs", grup_B, "sonuçları çekiliyor... (kriter:", kriter, ")\n")
  res <- results(dds, contrast = c("condition", grup_A, grup_B), alpha = p_siniri)
  res_tablo <- as.data.frame(res)
  res_tablo <- res_tablo[!is.na(res_tablo[[kriter]]), ]
  
  if (is.null(lfc_siniri)) {
    anlamli <- res_tablo[res_tablo[[kriter]] < p_siniri, ]
  } else {
    anlamli <- res_tablo[res_tablo[[kriter]] < p_siniri & abs(res_tablo$log2FoldChange) > lfc_siniri, ]
  }
  
  # Ensembl ID'lerindeki versiyon eklerini temizliyoruz
  anlamli$Gene_ID <- sub("\\.[0-9]+$", "", rownames(anlamli))
  anlamli$Symbol  <- sembol_getir(anlamli$Gene_ID)
  
  anlamli <- anlamli[, c("Gene_ID", "Symbol", "log2FoldChange", "pvalue", "padj")]
  anlamli <- anlamli[order(anlamli[[kriter]]), ]
  dosya_adi <- paste0("vGluT2_", grup_A, "_vs_", grup_B, "_DESeq2_RawData_", kriter, ".csv")
  write.csv(anlamli, dosya_adi, row.names = FALSE)
  
  cat("Kaydedilen dosya:", dosya_adi, "\n")
  cat("Anlamlı gen sayısı (", kriter, "<", p_siniri, "):", nrow(anlamli), "\n")
  if (nrow(anlamli) == 0) {
    cat("  UYARI: Bu eşikte hiç anlamlı gen bulunamadı.\n")
  } else {
    cat("  Yukarı regüle (log2FC > 0):", sum(anlamli$log2FoldChange > 0), "\n")
    cat("  Aşağı regüle (log2FC < 0):", sum(anlamli$log2FoldChange < 0), "\n")
  }
  return(anlamli)
}



#Üç Karşılaştırmanın Yapılması
cat("\nGruplar arası karşılaştırmalar başlıyor...\n")
cjd_control_genleri <- sonuc_cikar_ve_kaydet("CJD", "Control", kriter = "pvalue")
ffi_control_genleri <- sonuc_cikar_ve_kaydet("FFI", "Control", kriter = "pvalue")
cjd_ffi_genleri     <- sonuc_cikar_ve_kaydet("CJD", "FFI", kriter = "pvalue")



#CJD ve FFI Yanıtlarının Örtüşmesi (Kontrol'e Göre)
cjd_ids <- cjd_control_genleri$Gene_ID
ffi_ids <- ffi_control_genleri$Gene_ID

ortak      <- intersect(cjd_ids, ffi_ids)   # paylaşılan prion cevabı
sadece_cjd <- setdiff(cjd_ids, ffi_ids)     # CJD'ye özel
sadece_ffi <- setdiff(ffi_ids, cjd_ids)     # FFI'ye özel

cat("\nÖrtüşme -> Ortak:", length(ortak),
    "| Sadece CJD:", length(sadece_cjd),
    "| Sadece FFI:", length(sadece_ffi), "\n")

writeLines(sadece_cjd, "vGluT2_sadece_CJD_DESeq2.txt")
writeLines(sadece_ffi, "vGluT2_sadece_FFI_DESeq2.txt")
writeLines(ortak,      "vGluT2_ortak_DESeq2.txt")

cat("\nTüm analizler (DESeq2 + raw count) tamamlandı!\n")
cat("Üretilen dosyalar:\n")
cat(" - GSE198063_vGluT2_before_after_normalizasyon.png\n")
cat(" - GSE198063_vGluT2_PCA_normalizasyon_sonrasi.png\n")
cat(" - vGluT2_*_DESeq2_RawData_*.csv\n")
cat(" - vGluT2_ortak_DESeq2.txt / vGluT2_sadece_CJD_DESeq2.txt / vGluT2_sadece_FFI_DESeq2.txt\n")
