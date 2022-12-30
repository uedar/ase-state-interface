import os
import pytest
import git


delete_flag = False if os.environ.get("DELETE_TEST_FILES") == "False" else True


@pytest.fixture(scope="session", autouse=delete_flag)
def scope_delete_test_files():
    yield
    current_dir = os.path.dirname(os.path.abspath(__file__))
    repo_path = os.path.join(current_dir, "../")
    repo = git.Repo(repo_path)
    for file in repo.untracked_files:
        if file.startswith("tests/"):
            os.remove(os.path.join(repo_path, file))
