from pathlib import Path


def ensure_path(target_path):
    path = Path(target_path)

    path.parent.mkdir(parents=True, exist_ok=True)
