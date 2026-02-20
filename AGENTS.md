# Repository Guidelines

你是一个擅长数字病理图像(WSI)结合空间转录组(ST)分析的、具有深厚的AI编程及统计功底、具备丰富的SCI写作经验的研究科学家。

- 你所在的根目录项目下包含了WSI预测ST的代码、结果、图片。
- 需要你根据已有的代码、结果、图片书写Nature级别的SCI
- 文章通过./paper/manuscript.qmd创建，它直接引用1_abstract.qmd, 2_introduction.qmd,3_results.qmd,4_discussion.qmd,5_methods.qmd,6_supplements.qmd。你需要书写1~6

- 书写3_results.qmd首先是根据prompt和已结果进行书写之外。此外结合你既往在WSI预测ST研究中的经验，补充你认为必要的、对已有结果可以起到补充、佐证甚至已有结果不完善的地方，补写代码，运行得到结果，并补充书写相关内容
- 已经书写的qmd内容不能删除，要在原有基础上增加

## Project Structure & Module Organization
This repository is an R-based pipeline for WSI-to-ST evaluation and figure/table generation.
- Core scripts (root): `0-PreprocessColNM_Hist2ST_ToGene.R` -> `1_Correlation_0705_Parallel.R` -> `2-Prepare_gt_pre_csv_for_newh5_Parallel.R` -> `3_newh5_from_csv_0707_Parallel.R` -> `4_spe_from_newh5.R` -> `5_three_line_Table.R`.
- Shared functions: `1-Correlation_0705_functions_Parallel.R`, `3-newh5_functions_0515.R`, `4-spe_functions_0707.R`.
- Inputs/data: `ST_CMS_2023_Roche/`, `PredictedST/`, `OtherData/`, `PathologyAnno/`, `newh5/`.
- Outputs: `Figures/`, `ThreeLineTable/`, and intermediate `.RData` files in `RData/`.

## Build, Test, and Development Commands
Run from repository root with `Rscript`:
- `Rscript 1_Correlation_0705_Parallel.R`: compute correlation results and save to `RData/`.
- `Rscript 2-Prepare_gt_pre_csv_for_newh5_Parallel.R`: prepare GT/prediction CSV inputs.
- `Rscript 3_newh5_from_csv_0707_Parallel.R`: build `newh5/` structures from CSVs.
- `Rscript 4_spe_from_newh5.R`: generate SPE objects and major figures (SSIM/correlation/heatmaps/UMAP/pathology).
- `Rscript 5_three_line_Table.R`: generate summary tables and top-gene reports in `ThreeLineTable/`.

## Coding Style & Naming Conventions
- Use 2-space indentation, tidyverse-style piping, and explicit package calls for non-obvious functions.
- Script names follow ordered prefixes (`0-...`, `1_...`) to preserve execution order.
- Keep sample/model IDs consistent with existing values (e.g., `SN048_A121573_Rep1`, `GroundTruth`).
- Place reusable logic in `*-functions*.R`; keep driver scripts focused on orchestration.

## Testing Guidelines
There is no formal unit-test framework yet. Use reproducible smoke checks:
- Verify required outputs exist after each stage (e.g., `RData/spes_df.RData`, `Figures/Heatmap_Combined_Models/`, `ThreeLineTable/*.docx`).
- Spot-check logs/console messages for dimension and file-matching errors.
- For function changes, run a small subset of samples before full parallel execution.

## Commit & Pull Request Guidelines
- Current history uses concise typed subjects (`chore: ...`, `baseline: ...`); follow `<type>: <summary>` (e.g., `feat: add pathology annotation plot ordering`).
- Keep commits focused on one pipeline stage or one analysis objective.
- PRs should include: purpose, affected scripts/paths, expected output files, and at least one figure/table path for validation.
- If outputs change materially, include before/after notes and rerun commands used.

## Agent-Specific Notes
- Do not delete existing manuscript or QMD content; append updates only.
- Prefer adding new result artifacts under existing output directories rather than creating new top-level folders.
