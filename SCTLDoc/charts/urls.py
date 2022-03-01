from django.urls import path
from .views import HomePageView, AboutPageView, AboutChartsView

urlpatterns = [
    path("charts/", AboutChartsView.as_view(), name="charts"),
    path("about/", AboutPageView.as_view(), name="about"),
    path("", HomePageView.as_view(), name="home"),
]
