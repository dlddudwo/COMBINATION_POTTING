#include "pch.h"
#include "CombinationFeatureDual_CV_InputProvider.h"

SCombinationDualCvInput CCombinationFeatureDualCvInputProvider::GetInput(const std::string& a_panel_id, int a_cam_index, int a_pattern_index) const
{
	UNREFERENCED_PARAMETER(a_panel_id);
	UNREFERENCED_PARAMETER(a_cam_index);
	UNREFERENCED_PARAMETER(a_pattern_index);

	SCombinationDualCvInput input;
	input.valid = false;
	return input;
}
