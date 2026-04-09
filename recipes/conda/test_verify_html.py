"""Verify HTML report contains expected content."""
with open("test_report.html") as f:
    content = f.read()
    assert "<!DOCTYPE html>" in content or "<html" in content
    assert "covsnap" in content.lower()
    assert "BRCA1" in content
    assert "hg38" in content
    print("All HTML report checks passed.")
