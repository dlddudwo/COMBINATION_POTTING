#pragma once

#include "FuncAlgorithm.h"

// Step 1: Structure-only refactor target.
// NOTE: This class is intentionally separated from CCombinationFeatureDual.
class CCombinationFeatureDual_CV
{
public:
	virtual Json::Value Run(Json::Value a_recipe, Json::Value& a_result_json);

private:
	struct SPatternTypeIndex
	{
		int omit = 0;
		int domit = 0;
		int black_domit = 0;
	};

	struct SRunContext
	{
		int cam_num = 0;
		int cam_index = 0;
		int pattern_count = 0;
		std::string panel_id;
	};

	SRunContext BuildRunContext(const Json::Value& a_recipe, const Json::Value& a_result_json) const;
	SPatternTypeIndex FindPatternTypeIndex(const Json::Value& a_recipe, int a_cam_index) const;
	void ProcessPatterns(const Json::Value& a_recipe, Json::Value& a_result_json, const SRunContext& a_ctx, const SPatternTypeIndex& a_pattern_index);
};
