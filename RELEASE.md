# Release Process

This document describes how to release a new version of `cij`.

## Single Source of Truth

The version is maintained in **one place only**: `cij/version.py`

The `pyproject.toml` dynamically reads the version from this file at build time.

## Quick Release (Automated)

### Using bump-my-version

```bash
# Install bump-my-version (one-time setup)
pip install bump-my-version

# Bump version (automatically commits and tags)
bump-my-version bump patch  # 1.0.0-b4 -> 1.0.0-b5
bump-my-version bump minor  # 1.0.0-b4 -> 1.1.0
bump-my-version bump major  # 1.0.0-b4 -> 2.0.0
```

This automatically:
1. Updates `cij/version.py`
2. Creates a git commit with message "Bump version to X.Y.Z"
3. Creates a git tag `vX.Y.Z`

### Push and Release

```bash
# Push commits and tags
git push && git push --tags
```

Once the tag is pushed, GitHub Actions automatically:
1. Builds the package
2. Publishes to PyPI
3. Creates a GitHub Release with changelog

## Manual Release (Without bump-my-version)

If you prefer to bump versions manually:

### 1. Update version

Edit `cij/version.py`:
```python
__version__ = "1.0.1"  # Your new version
```

### 2. Commit and tag

```bash
git add cij/version.py
git commit -m "Bump version to 1.0.1"
git tag v1.0.1
```

### 3. Push to GitHub

```bash
git push && git push --tags
```

### 4. GitHub Actions takes over

The `.github/workflows/release.yml` workflow will automatically:
- Build the package with `uv`
- Upload to PyPI
- Create a GitHub Release with auto-generated changelog

## Pre-releases

For beta/alpha releases, use version strings like:
- `1.0.0-b1` (beta 1)
- `1.0.0-alpha1` (alpha 1)
- `1.0.0-rc1` (release candidate 1)

The GitHub Actions workflow will automatically mark these as pre-releases.

## Verification

After release, verify:

1. **PyPI**: Check https://pypi.org/project/cij/
2. **GitHub Releases**: Check https://github.com/MineralsCloud/cij/releases
3. **Installation**:
   ```bash
   pip install --upgrade cij
   python -c "import cij; print(cij.__version__)"
   ```

## Troubleshooting

### Build fails locally

```bash
# Test the build locally
python -m pip install build
python -m build
```

### Version not updating

Make sure `cij/version.py` is committed before tagging:
```bash
git status  # Should show nothing or only untracked files
```

### PyPI upload fails

Check that secrets are configured in GitHub:
- `PYPI_USERNAME`
- `PYPI_PASSWORD` (or `PYPI_API_TOKEN`)

## Configuration Files

- `cij/version.py` - Single source of truth for version
- `pyproject.toml` - Project metadata, dynamic version, and bump-my-version config
- `.github/workflows/release.yml` - Automated release workflow (uses uv)
