#!/bin/bash
set -euo pipefail

# 使用方法を表示する関数
usage() {
    echo "使用方法: $0 [サンプル名] [リードディレクトリ] [リファレンスゲノム] [出力プレフィックス]"
    echo ""
    echo "引数:"
    echo "  サンプル名         - 処理するサンプルの名前 (デフォルト: sample)"
    echo "                      このパターンを使って *_1.fastq.gz と *_2.fastq.gz を検索します"
    echo "  リードディレクトリ - FASTQファイルが格納されているディレクトリ (デフォルト: カレントディレクトリ)"
    echo "  リファレンスゲノム - 使用するリファレンスゲノムのパス (デフォルト: hg38.fa)"
    echo "  出力プレフィックス - 出力ファイルのプレフィックス (デフォルト: サンプル名と同じ)"
    echo ""
    echo "例:"
    echo "  $0                           # デフォルト値を使用"
    echo "  $0 MY_SAMPLE                 # サンプル名を指定"
    echo "  $0 MY_SAMPLE /path/to/reads  # サンプル名とディレクトリを指定"
    exit 1
}

# ヘルプオプションのチェック
if [[ "$#" -gt 0 && ("$1" == "-h" || "$1" == "--help") ]]; then
    usage
fi

# 入力ファイル設定
SAMPLE_NAME=${1:-"sample"}
READ_DIR=${2:-"."}
READ_PATTERN="${SAMPLE_NAME}*"

# リードファイルの自動検出
READ1=$(find ${READ_DIR} -name "${READ_PATTERN}*_1.fastq.gz" | head -n 1)
READ2=$(find ${READ_DIR} -name "${READ_PATTERN}*_2.fastq.gz" | head -n 1)

if [ -z "$READ1" ] || [ -z "$READ2" ]; then
    echo "エラー: リードファイルが見つかりません。パターン: ${READ_PATTERN}"
    usage
fi

REF=${3:-"hg38.fa"}
OUTPREFIX=${4:-$SAMPLE_NAME}
THREADS=${THREADS:-4}
JAVA_OPTS=${JAVA_OPTS:-"-Xmx4g"}
PICARD_CMD="picard"

echo "==== NGSリードマッピングパイプライン開始 ===="
echo "サンプル名: $SAMPLE_NAME"
echo "読み込むリード:"
echo "- R1: $READ1"
echo "- R2: $READ2"
echo "リファレンス: $REF"
echo "出力プレフィックス: $OUTPREFIX"

# 1. fastpでトリミング
echo "1. トリミング実行中..."
fastp -i ${READ1} -I ${READ2} \
      -o ${OUTPREFIX}_trimmed_R1.fastq.gz -O ${OUTPREFIX}_trimmed_R2.fastq.gz \
      --detect_adapter_for_pe \
      --thread ${THREADS} \
      -h ${OUTPREFIX}_fastp.html -j ${OUTPREFIX}_fastp.json

# 2. BWAインデックス確認
echo "2. BWAインデックス確認中..."
if [ ! -f "${REF}.bwt" ]; then
    echo "BWAインデックスが見つかりません。作成します..."
    bwa index ${REF}
else
    echo "インデックスあり。スキップ。"
fi

# 3. マッピング（Read Group付き）
echo "3. BWAマッピング中..."
bwa mem -t ${THREADS} \
  -R "@RG\tID:${OUTPREFIX}\tSM:${OUTPREFIX}\tPL:ILLUMINA" \
  ${REF} \
  ${OUTPREFIX}_trimmed_R1.fastq.gz ${OUTPREFIX}_trimmed_R2.fastq.gz > ${OUTPREFIX}.sam

# 4. SAM→BAM→ソート
echo "4. SAMをBAMへ変換＆ソート中..."
samtools view -@ ${THREADS} -bS ${OUTPREFIX}.sam | \
  samtools sort -@ ${THREADS} -o ${OUTPREFIX}.sorted.bam -

# 5. 重複マーク
echo "5. 重複リードをマーク中..."
${PICARD_CMD} MarkDuplicates \
    INPUT=${OUTPREFIX}.sorted.bam \
    OUTPUT=${OUTPREFIX}.dedup.bam \
    METRICS_FILE=${OUTPREFIX}.dedup_metrics.txt \
    REMOVE_DUPLICATES=false \
    ASSUME_SORTED=true \
    VALIDATION_STRINGENCY=LENIENT

# 6. インデックス作成
echo "6. インデックス作成中..."
samtools index ${OUTPREFIX}.dedup.bam

# 7. 統計出力
echo "7. BAM統計を出力中..."
samtools flagstat ${OUTPREFIX}.dedup.bam > ${OUTPREFIX}.flagstat.txt
samtools idxstats ${OUTPREFIX}.dedup.bam > ${OUTPREFIX}.idxstats.txt
samtools stats ${OUTPREFIX}.dedup.bam > ${OUTPREFIX}.stats.txt

# 終了メッセージ
echo "==== パイプライン完了 ===="
echo "出力ファイル:"
echo "- ${OUTPREFIX}.dedup.bam (重複マーク済み)"
echo "- ${OUTPREFIX}.dedup.bam.bai"
echo "- 統計: ${OUTPREFIX}.flagstat.txt / ${OUTPREFIX}.idxstats.txt / ${OUTPREFIX}.stats.txt"
echo "- fastp レポート: ${OUTPREFIX}_fastp.html / ${OUTPREFIX}_fastp.json"
echo "- Picard メトリクス: ${OUTPREFIX}.dedup_metrics.txt"
