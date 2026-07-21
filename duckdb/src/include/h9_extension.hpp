#pragma once

#include "duckdb.hpp"

namespace duckdb {

class H9Extension : public Extension {
public:
	void Load(ExtensionLoader &loader) override;
	std::string Name() override;
	std::string Version() const override;
};

} // namespace duckdb
