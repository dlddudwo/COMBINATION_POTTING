# CombinationFeatureDual_CV Golden Compare Test Draft

## Goal
Validate that `CCombinationFeatureDual_CV::Run` produces equivalent result-json outputs to `CCombinationFeatureDual::Run` before HALCON->OpenCV migration.

## Draft Procedure
1. Prepare same `recipe` + initialized `result_json` fixture for both classes.
2. Execute legacy path:
   - `CCombinationFeatureDual legacy;`
   - `Json::Value legacy_out = legacy.Run(recipe, legacy_result);`
3. Execute CV refactor path:
   - `CCombinationFeatureDual_CV cv;`
   - `Json::Value cv_out = cv.Run(recipe, cv_result);`
4. Compare recursively:
   - key existence (missing/extra)
   - numeric tolerance (epsilon: `1e-6`)
   - string/int exact equality
5. Emit per-PTN/per-DEFECT diff report.

## Priority Feature Keys
- `OmitScore`, `DOmitScore`, `BlackDOmitScore`
- `Omit_Ratio_64`, `DOmit_Ratio_64`, `Black_DOmit_Ratio_64`
- `OmitAvg_Pre`, `OmitAvg_Ori`, `DomitAvg_Pre`, `DomitAvg_Ori`
- `GrayAVG_Pre64`, `GrayMin_Pre64`, `GrayMax_Pre64`, `GraySTDEV_Pre64`
- `GrayOver135_Pre64` ~ `GrayOver170_Pre64`
- `GrayInner135_Pre` ~ `GrayInner240_Pre`

## Pass Criteria
- No missing keys.
- No extra keys.
- Numeric mismatches within tolerance only.
