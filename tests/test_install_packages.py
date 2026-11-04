"""Package source classification for the R setup script."""

from pathlib import Path
import re
import unittest


INSTALLER = Path(__file__).resolve().parents[1] / "pipeline" / "install_packages.R"


def package_vector(script, name):
    match = re.search(rf"{name}\s*<-\s*c\((.*?)\)", script, re.DOTALL)
    if match is None:
        raise AssertionError(f"Could not find {name} in {INSTALLER}")
    return re.findall(r'"([^"]+)"', match.group(1))


class InstallerPackageSourceTests(unittest.TestCase):
    def test_decontam_is_installed_from_bioconductor(self):
        script = INSTALLER.read_text(encoding="utf-8")

        self.assertNotIn("decontam", package_vector(script, "cran_pkgs"))
        self.assertIn("decontam", package_vector(script, "bioc_pkgs"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
