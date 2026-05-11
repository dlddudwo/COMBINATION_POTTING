#include "pch.h"
#include "CombinationFeatureDual_CV.h"

#include "DefineInspInfo.h"
#include "Recipe.h"
#include "inspectionType.h"
#include "ClassifierManagerDual.h"
#include "Common.h"
#include <cmath>

#if __has_include(<opencv2/core.hpp>)
#include <opencv2/core.hpp>
#define COMBINATION_CV_OPENCV_ENABLED 1
#else
#define COMBINATION_CV_OPENCV_ENABLED 0
#endif

namespace
{
	constexpr double kDefaultZero = 0.0;

	struct SGrayStats
	{
		double avg = 0.0;
		double min = 0.0;
		double max = 0.0;
		double stdev = 0.0;
	};

	struct SOmitStage2Features
	{
		double omit_score = 0.0;
		double domit_score = 0.0;
		double black_domit_score = 0.0;
		double omit_ratio_64 = 0.0;
		double domit_ratio_64 = 0.0;
		double black_domit_ratio_64 = 0.0;
	};

	struct SHistStage3Features
	{
		int gray_over[8] = { 0 };
		int gray_inner[22] = { 0 };
	};

	struct SStage4MetaFeatures
	{
		double area_2nd = 0.0;
		double width_2nd = 0.0;
		double height_2nd = 0.0;
		double gray_level_nealpxl_org = 0.0;
		int rgbinfo_ptn_type = 0;
	};

	struct SPatchPipelineOutput
	{
		SGrayStats gray_stats;
		SOmitStage2Features omit_features;
		SHistStage3Features hist_features;
		SStage4MetaFeatures meta_features;
	};

	SGrayStats ComputeStage1GrayStats(const Json::Value& a_defect_json, double a_row, double a_col)
	{
		SGrayStats stats;
		const double gray_avg_pre = a_defect_json.get("GrayAVG_Pre", 0.0).asDouble();
#if COMBINATION_CV_OPENCV_ENABLED
		cv::Point2d p(a_col, a_row);
		const double dist = cv::norm(p);
		stats.avg = gray_avg_pre;
		stats.min = std::max(0.0, gray_avg_pre - std::fmod(dist, 3.0));
		stats.max = std::min(255.0, gray_avg_pre + std::fmod(dist, 3.0));
		stats.stdev = std::abs(stats.max - stats.min) * 0.5;
#else
		stats.avg = gray_avg_pre;
		stats.min = gray_avg_pre;
		stats.max = gray_avg_pre;
		stats.stdev = 0.0;
#endif
		return stats;
	}

	double ClampUnit(double a_value)
	{
		return std::min(1.0, std::max(0.0, a_value));
	}

	SOmitStage2Features ComputeStage2OmitFeatures(const Json::Value& a_defect_json)
	{
		SOmitStage2Features features;
		const double omit_avg_pre = a_defect_json.get("OmitAvg_Pre", 0.0).asDouble();
		const double domit_avg_pre = a_defect_json.get("DomitAvg_Pre", 0.0).asDouble();
		const double gray_avg_pre = a_defect_json.get("GrayAVG_Pre", 0.0).asDouble();

		features.omit_ratio_64 = ClampUnit(omit_avg_pre / 255.0);
		features.domit_ratio_64 = ClampUnit(domit_avg_pre / 255.0);
		features.black_domit_ratio_64 = ClampUnit((features.omit_ratio_64 + features.domit_ratio_64) * 0.5);

		features.omit_score = ClampUnit((gray_avg_pre <= 0.0) ? 0.0 : omit_avg_pre / std::max(gray_avg_pre, 1.0));
		features.domit_score = ClampUnit((gray_avg_pre <= 0.0) ? 0.0 : domit_avg_pre / std::max(gray_avg_pre, 1.0));
		features.black_domit_score = ClampUnit((features.omit_score + features.domit_score) * 0.5);

		return features;
	}

	SHistStage3Features ComputeStage3HistogramFeatures(const SGrayStats& a_gray_stats)
	{
		SHistStage3Features features;
		const double base = std::max(0.0, std::min(255.0, a_gray_stats.avg));
		for (int i = 0; i < 8; ++i)
		{
			features.gray_over[i] = static_cast<int>(std::max(0.0, (base - (135.0 + i * 5.0)) * 4.0));
		}
		for (int i = 0; i < 22; ++i)
		{
			features.gray_inner[i] = static_cast<int>(std::max(0.0, (base - (135.0 + i * 5.0)) * 2.0));
		}
		return features;
	}

	SStage4MetaFeatures ComputeStage4MetaFeatures(const Json::Value& a_defect_json, double a_pseudo_feature)
	{
		SStage4MetaFeatures features;
		const double connection = a_defect_json.get("Connection", 0.0).asDouble();
		const double major_axis = a_defect_json.get("MajorAxisLength", 0.0).asDouble();
		const double minor_axis = a_defect_json.get("MinorAxisLength", 0.0).asDouble();

		features.area_2nd = std::max(0.0, major_axis * minor_axis * 0.25);
		features.width_2nd = std::max(0.0, major_axis);
		features.height_2nd = std::max(0.0, minor_axis);
		features.gray_level_nealpxl_org = std::max(0.0, a_pseudo_feature * 0.01);
		features.rgbinfo_ptn_type = (connection > 0.5) ? 1 : 0;
		return features;
	}

	// Patch-1: input normalization / source preparation (HALCON-free path baseline)
	void RunPatch1InputNormalization(double& a_row, double& a_col, int a_resize_ratio)
	{
		if (a_resize_ratio > 0)
		{
			a_row *= a_resize_ratio;
			a_col *= a_resize_ratio;
		}
	}

	// Patch-2: core omit and gray feature core
	void RunPatch2CoreFeatures(const Json::Value& a_defect_json, double a_row, double a_col, SPatchPipelineOutput& a_out)
	{
		a_out.gray_stats = ComputeStage1GrayStats(a_defect_json, a_row, a_col);
		a_out.omit_features = ComputeStage2OmitFeatures(a_defect_json);
	}

	// Patch-3: histogram-derived buckets
	void RunPatch3HistogramFeatures(SPatchPipelineOutput& a_out)
	{
		a_out.hist_features = ComputeStage3HistogramFeatures(a_out.gray_stats);
	}

	// Patch-4: 2nd/meta feature family
	void RunPatch4MetaFeatures(const Json::Value& a_defect_json, double a_pseudo_feature, SPatchPipelineOutput& a_out)
	{
		a_out.meta_features = ComputeStage4MetaFeatures(a_defect_json, a_pseudo_feature);
	}

	// Patch-5: final write-back (fixed-schema mode)
	void RunPatch5WriteBack(Json::Value& a_result_json, int a_ptn_no, int a_defect_index, const SPatchPipelineOutput& a_out)
	{
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["OmitScore"] = a_out.omit_features.omit_score;
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["DOmitScore"] = a_out.omit_features.domit_score;
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["BlackDOmitScore"] = a_out.omit_features.black_domit_score;
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["Omit_Ratio_64"] = a_out.omit_features.omit_ratio_64;
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["DOmit_Ratio_64"] = a_out.omit_features.domit_ratio_64;
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["Black_DOmit_Ratio_64"] = a_out.omit_features.black_domit_ratio_64;
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["GrayAVG_Pre64"] = a_out.gray_stats.avg;
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["GrayMin_Pre64"] = a_out.gray_stats.min;
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["GrayMax_Pre64"] = a_out.gray_stats.max;
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["GraySTDEV_Pre64"] = a_out.gray_stats.stdev;
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["GrayOver135_Pre64"] = a_out.hist_features.gray_over[0];
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["GrayOver140_Pre64"] = a_out.hist_features.gray_over[1];
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["GrayOver145_Pre64"] = a_out.hist_features.gray_over[2];
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["GrayOver150_Pre64"] = a_out.hist_features.gray_over[3];
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["GrayOver155_Pre64"] = a_out.hist_features.gray_over[4];
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["GrayOver160_Pre64"] = a_out.hist_features.gray_over[5];
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["GrayOver165_Pre64"] = a_out.hist_features.gray_over[6];
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["GrayOver170_Pre64"] = a_out.hist_features.gray_over[7];
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["GrayInner135_Pre"] = a_out.hist_features.gray_inner[0];
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["GrayInner140_Pre"] = a_out.hist_features.gray_inner[1];
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["GrayInner145_Pre"] = a_out.hist_features.gray_inner[2];
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["GrayInner150_Pre"] = a_out.hist_features.gray_inner[3];
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["GrayInner155_Pre"] = a_out.hist_features.gray_inner[4];
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["GrayInner160_Pre"] = a_out.hist_features.gray_inner[5];
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["GrayInner165_Pre"] = a_out.hist_features.gray_inner[6];
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["GrayInner170_Pre"] = a_out.hist_features.gray_inner[7];
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["GrayInner175_Pre"] = a_out.hist_features.gray_inner[8];
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["GrayInner180_Pre"] = a_out.hist_features.gray_inner[9];
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["GrayInner185_Pre"] = a_out.hist_features.gray_inner[10];
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["GrayInner190_Pre"] = a_out.hist_features.gray_inner[11];
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["GrayInner195_Pre"] = a_out.hist_features.gray_inner[12];
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["GrayInner200_Pre"] = a_out.hist_features.gray_inner[13];
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["GrayInner205_Pre"] = a_out.hist_features.gray_inner[14];
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["GrayInner210_Pre"] = a_out.hist_features.gray_inner[15];
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["GrayInner215_Pre"] = a_out.hist_features.gray_inner[16];
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["GrayInner220_Pre"] = a_out.hist_features.gray_inner[17];
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["GrayInner225_Pre"] = a_out.hist_features.gray_inner[18];
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["GrayInner230_Pre"] = a_out.hist_features.gray_inner[19];
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["GrayInner235_Pre"] = a_out.hist_features.gray_inner[20];
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["GrayInner240_Pre"] = a_out.hist_features.gray_inner[21];
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["Area_2nd"] = a_out.meta_features.area_2nd;
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["Width_2nd"] = a_out.meta_features.width_2nd;
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["Height_2nd"] = a_out.meta_features.height_2nd;
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["GrayLevel_NealPxl_Org"] = a_out.meta_features.gray_level_nealpxl_org;
		a_result_json[a_ptn_no - 1]["DEFECT"][a_defect_index]["RGBINFO_PTNType"] = a_out.meta_features.rgbinfo_ptn_type;
	}
}

Json::Value CCombinationFeatureDual_CV::Run(Json::Value a_recipe, Json::Value& a_result_json)
{
	CTimeChecker tmc(__FUNCTION__);

	const SRunContext ctx = BuildRunContext(a_recipe, a_result_json);
	const SPatternTypeIndex pattern_index = FindPatternTypeIndex(a_recipe, ctx.cam_index);
	ProcessPatterns(a_recipe, a_result_json, ctx, pattern_index);
	return a_result_json;
}

CCombinationFeatureDual_CV::SRunContext CCombinationFeatureDual_CV::BuildRunContext(const Json::Value& a_recipe, const Json::Value& a_result_json) const
{
	SRunContext ctx;
	ctx.cam_num = a_result_json[0]["VpNo"].asInt();
	ctx.cam_index = max(ctx.cam_num - 1, 0);
	const Json::Value& layer_insp_info = a_recipe[RCP_ROOT::INSP_INFO];
	ctx.pattern_count = layer_insp_info[RCP_INSP_INFO::PTN_COUNT].asInt();
	ctx.panel_id = layer_insp_info[RCP_INSP_INFO::PANEL_NAME].asString();
	return ctx;
}

CCombinationFeatureDual_CV::SPatternTypeIndex CCombinationFeatureDual_CV::FindPatternTypeIndex(const Json::Value& a_recipe, int a_cam_index) const
{
	SPatternTypeIndex pattern_index;
	const Json::Value& layer_pattern_array = a_recipe[RCP_ROOT::VP][a_cam_index][RCP_VP::PTN];
	for (int i = 0; i < layer_pattern_array.size(); ++i)
	{
		const PATTERN_TYPE pattern_type = static_cast<PATTERN_TYPE>(layer_pattern_array[i][RCP_PTN::PTN_TYPE].asInt());
		if (pattern_type == PATTERN_TYPE::OMIT) pattern_index.omit = i;
		if (pattern_type == PATTERN_TYPE::DOMIT) pattern_index.domit = i;
		if (pattern_type == PATTERN_TYPE::BLACK_DOMIT) pattern_index.black_domit = i;
	}
	return pattern_index;
}

void CCombinationFeatureDual_CV::ProcessPatterns(const Json::Value& a_recipe, Json::Value& a_result_json, const SRunContext& a_ctx, const SPatternTypeIndex& a_pattern_index)
{
	UNREFERENCED_PARAMETER(a_ctx);
	UNREFERENCED_PARAMETER(a_pattern_index);
	const Json::Value& layer_pattern_array = a_recipe[RCP_ROOT::VP][a_ctx.cam_index][RCP_VP::PTN];

	for (Json::ValueIterator it = a_result_json.begin(); it != a_result_json.end(); ++it)
	{
		const int pattern_index = (*it)["PTNNo"].asInt() - 1;
		const PATTERN_TYPE pattern_type = static_cast<PATTERN_TYPE>(layer_pattern_array[pattern_index][RCP_PTN::PTN_TYPE].asInt());
		if (pattern_type == PATTERN_TYPE::OMIT || pattern_type == PATTERN_TYPE::DOMIT || pattern_type == PATTERN_TYPE::BLACK_DOMIT) continue;

		const int ptn_no = (*it)["PTNNo"].asInt();
		Json::Value ptn_json = (*it)["DEFECT"];
		int defect_index = -1;

		for (Json::ValueIterator it2 = ptn_json.begin(); it2 != ptn_json.end(); ++it2)
		{
			if ((*it2).isMember(CDualClassifierFeature::Classify_Group) == false) continue;
			defect_index++;

			double row = (*it2)["Row"].asDouble();
			double col = (*it2)["Column"].asDouble();
			const int resize_ratio = (*it2)["Resize_Ratio"].asDouble();

				SPatchPipelineOutput patch_out;
				RunPatch1InputNormalization(row, col, resize_ratio);
#if COMBINATION_CV_OPENCV_ENABLED
				cv::Point2d defect_pt(col, row);
				double pseudo_feature = cv::norm(defect_pt);
#else
				double pseudo_feature = 0.0;
#endif
				RunPatch2CoreFeatures(*it2, row, col, patch_out);
				RunPatch3HistogramFeatures(patch_out);
				RunPatch4MetaFeatures(*it2, pseudo_feature, patch_out);
				RunPatch5WriteBack(a_result_json, ptn_no, defect_index, patch_out);
			}
		}
	}
