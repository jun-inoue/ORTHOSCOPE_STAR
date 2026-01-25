
#!/usr/bin/env bash
set -euo pipefail

# 使い方:
#   bash 000push.sh [FILE ...] [-m "commit message"]
# 例:
#   bash 000push.sh images/SpeciesTree_Urochordata.pdf
#   bash 000push.sh outdir_AA/*.txt -m "Update AA txts"
#   bash 000push.sh -m "Tweak README" README.md scripts/*.sh

# --- ヘルパー: 文字列の存在チェック ---
has_arg() { local x="$1"; shift; for a in "$@"; do [[ "$a" == "$x" ]] && return 0; done; return 1; }

# --- Gitリポジトリ確認 & ルート取得 ---
if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  echo "エラー: このディレクトリはGitリポジトリ内ではありません。" >&2
  exit 1
fi
REPO_ROOT="$(git rev-parse --show-toplevel)"

# --- 現在ブランチ & リモート検出 ---
CURRENT_BRANCH="$(git rev-parse --abbrev-ref HEAD)"
# もしdetached HEADなら既定ブランチを推定
if [[ "$CURRENT_BRANCH" == "HEAD" ]]; then
  if git show-ref --verify --quiet refs/heads/main; then
    CURRENT_BRANCH="main"
  elif git show-ref --verify --quiet refs/heads/master; then
    CURRENT_BRANCH="master"
  else
    echo "エラー: ブランチ名を特定できません（detached HEAD）。-m付きでブランチにcheckoutしてから実行してください。" >&2
    exit 1
  fi
fi

# 既定リモートはorigin。なければ最初のfetch可能リモートを選択
REMOTE="origin"
if ! git remote | grep -qx "$REMOTE"; then
  if git remote | grep -q .; then
    REMOTE="$(git remote | head -n1)"
    echo "注意: 'origin' が無いので '${REMOTE}' を使用します。"
  else
    echo "エラー: リモートが設定されていません。'git remote add origin ...' を先に実行してください。" >&2
    exit 1
  fi
fi

# --- 引数解析（-m メッセージ対応 / 複数ファイル可） ---
COMMIT_MSG=""
FILES=()

if has_arg "-m" "$@"; then
  # -m の後ろをコミットメッセージとして取り出す
  # シンプルなパーサ（複数 -m は非対応）
  SEEN_M=false
  for arg in "$@"; do
    if $SEEN_M; then
      COMMIT_MSG="$arg"
      SEEN_M=false
      continue
    fi
    if [[ "$arg" == "-m" ]]; then
      SEEN_M=true
      continue
    fi
    # -m とメッセージ以外はファイル候補に
    FILES+=("$arg")
  done
else
  # -m が無い場合は全引数をファイル候補
  for arg in "$@"; do
    FILES+=("$arg")
  done
fi

# --- デフォルト対象: README.md ---
if [[ "${#FILES[@]}" -eq 0 ]]; then
  FILES=("README.md")
fi

# --- glob展開（存在しないglobはそのまま文字列になるので補正） ---
EXPANDED_FILES=()
for f in "${FILES[@]}"; do
  shopt -s nullglob
  matches=($f)
  shopt -u nullglob
  if [[ "${#matches[@]}" -eq 0 ]]; then
    matches=("$f")
  fi
  EXPANDED_FILES+=("${matches[@]}")
done

# --- ファイル存在 & リポジトリ内パスチェック ---
STAGE_FILES=()
for f in "${EXPANDED_FILES[@]}"; do
  if [[ ! -e "$f" ]]; then
    echo "エラー: 指定されたファイルが見つかりません: $f" >&2
    exit 1
  fi
  # 絶対パス化してリポジトリ内か検証
  abs_f="$(cd "$(dirname "$f")" && pwd)/$(basename "$f")"
  case "$abs_f" in
    "$REPO_ROOT"/*) ;;
    *)
      echo "エラー: リポジトリ外のファイルは追加できません: $f" >&2
      exit 1
      ;;
  esac
  STAGE_FILES+=("$f")
done

# --- 追加 & コミット ---
git add -- "${STAGE_FILES[@]}"

if [[ -z "$COMMIT_MSG" ]]; then
  if [[ "${#STAGE_FILES[@]}" -eq 1 ]]; then
    base="$(basename "${STAGE_FILES[0]}")"
    COMMIT_MSG="Update ${base}"
  else
    COMMIT_MSG="Update ${#STAGE_FILES[@]} files"
  fi
fi

# 変更が無い場合はコミットをスキップ
if git diff --cached --quiet; then
  echo "ステージ済みの変更がありません（コミットはスキップされました）。"
else
  git commit -m "$COMMIT_MSG" -- "${STAGE_FILES[@]}"
fi

# --- push ---
git push "$REMOTE" "$CURRENT_BRANCH"
echo "✅ Push 完了: ${REMOTE} ${CURRENT_BRANCH}"
